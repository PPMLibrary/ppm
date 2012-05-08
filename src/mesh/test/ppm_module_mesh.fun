test_suite ppm_module_mesh

use ppm_module_mesh_typedef
use ppm_module_topo_typedef
use ppm_module_field_typedef
use ppm_module_mktopo

#ifdef __MPI
    INCLUDE "mpif.h"
#endif

integer, parameter              :: debug = 0
integer, parameter              :: mk = kind(1.0d0) !kind(1.0e0)
real(mk),parameter              :: pi = ACOS(-1._mk)
integer,parameter               :: ndim=3
integer                         :: decomp,assig,tolexp
integer                         :: info,comm,rank,nproc
real(mk)                        :: tol
integer                         :: topoid=-1
real(mk),dimension(:  ),pointer :: min_phys => NULL()
real(mk),dimension(:  ),pointer :: max_phys => NULL()

integer, dimension(:  ),pointer :: ighostsize => NULL()
real(mk)                        :: sca_ghostsize

integer                         :: i,j,k
integer                         :: nsublist
integer, dimension(:  ),pointer :: isublist => NULL()
integer, dimension(2*ndim)      :: bcdef
real(mk),dimension(:  ),pointer :: cost => NULL()
integer, dimension(:  ),pointer :: nm => NULL()
real(mk),dimension(:  ),pointer :: h => NULL()
type(ppm_t_topo),       pointer :: topo => NULL()

type(ppm_t_equi_mesh)            :: Mesh1,Mesh2
integer                          :: ipatch,isub
class(ppm_t_subpatch_),POINTER   :: p => NULL()

integer                          :: mypatchid
real(mk),dimension(2*ndim)       :: my_patch
real(mk),dimension(ndim)         :: offset

real(mk),dimension(:,:),pointer  :: field2d_1,field2d_2
real(mk),dimension(:,:,:),pointer:: field3d_1,field3d_2
real(mk),dimension(:,:,:,:),pointer:: field4d_1,field4d_2

!---------------- init -----------------------

    init

        use ppm_module_topo_typedef
        use ppm_module_init
        
        allocate(min_phys(ndim),max_phys(ndim),&
            &         ighostsize(ndim),nm(ndim),h(ndim))
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 1.0_mk
        ighostsize(1:ndim) = 2
        bcdef(1:2*ndim) = ppm_param_bcdef_periodic
        tolexp = -12

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

!----------------------------------------------

!---------------- finalzie --------------------


    finalize
        use ppm_module_finalize

        call ppm_finalize(info)

        deallocate(min_phys,max_phys,ighostsize,h,nm)

    end finalize

!----------------------------------------------

!------------- setup --------------------------

    setup

    end setup
!----------------------------------------------
        

!--------------- teardown ---------------------
    teardown
        NULLIFY(topo)
    end teardown
!----------------------------------------------

    test mesh_create_destroy
        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        assig  = ppm_param_assign_internal
        topoid = 0
        sca_ghostsize = 0.05_mk 
        call ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,    &
            &               bcdef,sca_ghostsize,cost,info)
        Assert_Equal(info,0)

        Nm = 125
        offset = 0._mk
        call Mesh1%create(topoid,offset,info,Nm=Nm)
        Assert_Equal(info,0)
        call Mesh1%destroy(info)
        Assert_Equal(info,0)

        h = (max_phys-min_phys)/Nm
        call Mesh1%create(topoid,offset,info,h=h)
        Assert_Equal(info,0)
        call Mesh1%destroy(info)
        Assert_Equal(info,0)
    end test

    test mesh_add_patches
        Nm = 129
        call Mesh1%create(topoid,offset,info,Nm=Nm)
        Assert_Equal(info,0)

        if (ndim .eq. 2) then

        mypatchid = 1
        my_patch(1:4) = (/0.5,0.1,5.1,10.0/)

        call Mesh1%def_patch(my_patch,info,mypatchid) 
        Assert_Equal(info,0)
        Assert_True(associated(Mesh1%subpatch))

        ipatch = 0
        topo => ppm_topo(Mesh1%topoid)%t
        do i = 1,topo%nsublist
            isub = topo%isublist(i)
            if (all(my_patch(1:ndim).LT.topo%max_subd(1:ndim,isub)) .AND. &
                all(my_patch(ndim+1:2*ppm_dim).GT.topo%min_subd(1:ndim,isub)))&
                then
                !count one subpatch
                ipatch = ipatch + 1
            endif
        enddo   

        isub = 0
        p => Mesh1%subpatch%begin()
        DO WHILE (ASSOCIATED(p))
            isub = isub+1
            p => Mesh1%subpatch%next()
        ENDDO
        Assert_Equal(isub,ipatch)

        mypatchid = 3
        my_patch(1:4) = (/0.1,0.0,0.3,0.8/)

        call Mesh1%def_patch(my_patch,info,mypatchid) 
        Assert_Equal(info,0)

        
        !find subpatches from patch 1
        DO i=1,Mesh1%patch%vec(1)%t%nsubpatch
            p => Mesh1%patch%vec(1)%t%subpatch(i)%t
            !write(*,*) 'subp no ',i,' for patch ',Mesh1%patch%vec(1)%t%patchid
        ENDDO
        p=>NULL()

        else
            write(*,*) 'TEST NOT RUNNING in 3D YET!!!'
        endif


        call Mesh1%destroy(info)
        Assert_Equal(info,0)

        Assert_False(associated(Mesh1%subpatch))
    end test

    test mesh_add_many_patches
        Nm = 129
        call Mesh1%create(topoid,offset,info,Nm=Nm)
        Assert_Equal(info,0)

        if (ndim .eq. 2) then

        mypatchid = 0
        DO i = 1,Nm(1)/4
            DO j = 1,Nm(2)/4
                mypatchid = mypatchid + 1
                my_patch(1:4) = (/h(1)*i,h(2)*j,h(1)*(i+4),h(2)*(j+4)/)
                call Mesh1%def_patch(my_patch,info,mypatchid) 
                Assert_Equal(info,0)
                Assert_True(associated(Mesh1%subpatch))
            ENDDO
        ENDDO

        else
            write(*,*) 'TEST NOT RUNNING in 3D YET!!!'
        endif
        
        call Mesh1%destroy(info)
        Assert_Equal(info,0)

        Assert_False(associated(Mesh1%subpatch))
    end test

    test mesh_add_patch_uniform
        !testing for a single patch that covers the whole domain
        Nm = 129
        Nm(ndim) = 65
        call Mesh1%create(topoid,offset,info,Nm=Nm)
            Assert_Equal(info,0)
        topo => ppm_topo(Mesh1%topoid)%t

        mypatchid = 1
        my_patch(1:ndim)        = min_phys(1:ndim)
        my_patch(ndim+1:2*ndim) = max_phys(1:ndim)

        write(*,*) 'adding patch : '
        write(*,*) my_patch(1:ndim)
        write(*,*) my_patch(ndim+1:2*ndim)

        call Mesh1%def_patch(my_patch,info,mypatchid) 
            Assert_Equal(info,0)

        Assert_True(associated(Mesh1%subpatch))

        isub = 0
        p => Mesh1%subpatch%begin()
        DO WHILE (ASSOCIATED(p))
            isub = isub+1
            p => Mesh1%subpatch%next()
        ENDDO
        Assert_Equal(isub,topo%nsublist)

        call Mesh1%destroy(info)
            Assert_Equal(info,0)

        !test the wrapper routine, which does the same thing
        call Mesh1%create(topoid,offset,info,Nm=Nm)
            Assert_Equal(info,0)
        call Mesh1%def_uniform(info)
            Assert_Equal(info,0)
        call Mesh1%destroy(info)
            Assert_Equal(info,0)
    end test

    test mesh_mappings_uniform
        !testing mappings for a single patch that covers the whole domain
        use ppm_module_io_vtk
        type(ppm_t_field) :: Field1
        real(ppm_kind_double),dimension(ndim) :: pos
        procedure(my_init_function), pointer :: init_f => NULL()

        init_f => my_init_function
        Nm = 129
        Nm(ndim) = 65
        Nm = 5
        call Mesh1%create(topoid,offset,info,Nm=Nm)
            Assert_Equal(info,0)
        call Mesh1%def_uniform(info)
            Assert_Equal(info,0)

        call Field1%create(3,info,name='testField',init_func=init_f) 
            Assert_Equal(info,0)
        call Field1%discretize_on(Mesh1,info)
            Assert_Equal(info,0)

        p => Mesh1%subpatch%begin()

        do while (ASSOCIATED(p))
            call p%get_field(field4d_1,Field1,info)
            do i = 1,p%nnodes(1)
                do j = 1,p%nnodes(2)
                    do k = 1,p%nnodes(3)
                        pos = p%get_pos(i,j)
                        field4d_1(1,i,j,k) = cos(2._mk*pi*pos(1))
                        field4d_1(2,i,j,k) = sin(2._mk*pi*pos(2))
                        field4d_1(3,i,j,k) = sin(2._mk*pi*pos(3))
                    enddo
                enddo
            enddo
            p => Mesh1%subpatch%next()
        enddo

        call ppm_vtk_fields('testvtk',Mesh1,info)
            Assert_Equal(info,0)

        !call Mesh1%map_ghost_get(info,ighostsize)
        !    Assert_Equal(info,0)

        call Mesh1%map_ghost_get(info,ighostsize)
            Assert_Equal(info,0)

        call Field1%map_ghost_push(Mesh1,info)
            Assert_Equal(info,0)
        write(*,*) 'PUSH OK'

        call Mesh1%map_send(info)
            Assert_Equal(info,0)
        write(*,*) 'SEND OK'

        call Field1%map_ghost_pop(Mesh1,info)
            Assert_Equal(info,0)
        write(*,*) 'POP OK'


        call Mesh1%destroy(info)
            Assert_Equal(info,0)
    end test

    FUNCTION my_init_function(x) RESULT(val)
        REAL(ppm_kind_double) :: val
        REAL(ppm_kind_double),DIMENSION(1:ndim),INTENT(IN) :: x

        val = sum(x)
    END FUNCTION
!============ Test cases ======================
!    test mesh_define
!        ! a simplistic test for checking if mesh_define is working
!
!        use ppm_module_topo_typedef
!        use ppm_module_mktopo
!        use ppm_module_map_field
!        use ppm_module_map_field_global
!        use ppm_module_mesh_define
!        use ppm_module_topo_get
!        use ppm_module_topo_check
!
!
!        allocate(nm(ndim),stat=info)
!        do i=1,ndim
!            nm(i) = 32*nproc
!        enddo
!
!        !----------------
!        ! make topology
!        !----------------
!        decomp = ppm_param_decomp_cuboid
!        !decomp = ppm_param_decomp_xpencil
!        assig  = ppm_param_assign_internal
!
!        topoid = 0
!        meshid = -1
!
!        call ppm_mktopo(topoid,meshid,xp,np,decomp,assig,min_phys,max_phys,    &
!        &               bcdef,ghostsize,cost,nm,info)
!
!        call ppm_topo_get_meshinfo(topoid,meshid,nm,istart,ndata,maxndata,&
!                        isublist,nsublist,info)
!        
!        allocate(field((1-ghostsize(1)):(maxndata(1)+ghostsize(1)),  &
!        &        (1-ghostsize(2)):(maxndata(2)+ghostsize(2)),nsublist),&
!        &        stat=info) ! 2d
!        allocate(field_ref((1-ghostsize(1)):(maxndata(1)+ghostsize(1)),  &
!        &        (1-ghostsize(2)):(maxndata(2)+ghostsize(2)),nsublist),&
!        &        stat=info) ! 2d
!        
!!        allocate(field((1-ghostsize(1)):(maxndata(1)+ghostsize(1)),  &
!!        &        (1-ghostsize(2)):(maxndata(2)+ghostsize(2)), &
!!        &        (1-ghostsize(3)):(maxndata(3)+ghostsize(3)),nsublist),&
!!        &       stat=info) ! 3d
!!        allocate(field_ref((1-ghostsize(1)):(maxndata(1)+ghostsize(1)),  &
!!        &        (1-ghostsize(2)):(maxndata(2)+ghostsize(2)), &
!!        &        (1-ghostsize(3)):(maxndata(3)+ghostsize(3)),nsublist),&
!!        &       stat=info) ! 3d
!        
!        do i=1,ndim
!            h(i) = (max_phys(i) - min_phys(i)) / real(ndata(i,1)-1,mk)
!        enddo
!       
!        meshid_ref = -1
!        nm_ref = nm
!        call ppm_mesh_define(topoid,meshid_ref,nm_ref,istart_ref,ndata_ref,info)
!        call ppm_map_field_global(topoid,topoid,meshid,meshid_ref,info)
!        call ppm_map_field_push(topoid,meshid,field,info)
!        call ppm_map_field_send(info)
!        call ppm_map_field_pop(topoid,meshid_ref,field_ref,ghostsize,info)
!
!
!        assert_equal(info,0)
!    end test


end test_suite
