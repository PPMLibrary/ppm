test_suite ppm_module_mesh



#ifdef __MPI
    INCLUDE "mpif.h"
#endif

integer, parameter              :: debug = 0
integer, parameter              :: mk = kind(1.0d0) !kind(1.0e0)
real(mk),parameter              :: pi = 3.1415926535897931_mk
integer,parameter               :: ndim=2
integer                         :: decomp,assig,tolexp
real(mk)                        :: tol
integer                         :: info,comm,rank,nproc
integer                         :: topoid
integer                         :: meshid,meshid_ref
real(mk),dimension(:  ),pointer :: min_phys => NULL()
real(mk),dimension(:  ),pointer :: max_phys => NULL()
real(mk),dimension(:  ),pointer :: h => NULL()
integer                         :: np
real(mk),dimension(:,:),pointer :: xp => NULL()
real(mk),dimension(:  ),pointer :: len_phys => NULL()
integer, dimension(:  ),pointer :: ghostsize => NULL()
integer                         :: i,j
integer                         :: nsublist
integer, dimension(:  ),pointer :: isublist => NULL()
integer, dimension(2*ndim)           :: bcdef
real(mk),dimension(:  ),pointer :: cost => NULL()
real(mk),dimension(:,:,:),pointer :: field => NULL()
real(mk),dimension(:,:,:),pointer :: field_ref => NULL()
integer, dimension(:  ),pointer :: nm => NULL()
integer, dimension(:  ),pointer :: nm_ref => NULL()
integer ,dimension(:,:),pointer :: istart => NULL()
integer ,dimension(:,:),pointer :: ndata => NULL()
integer ,dimension(:,:),pointer :: istart_ref => NULL()
integer ,dimension(:,:),pointer :: ndata_ref => NULL()
integer, dimension(ndim)        :: maxndata
integer                         :: seedsize
integer,  dimension(:),allocatable :: seed
real(mk), dimension(:),allocatable :: randnb
integer                          :: isymm = 0
logical                          :: lsymm = .false.,ok

!---------------- init -----------------------

    init

        use ppm_module_topo_typedef
        use ppm_module_init
        
        allocate(min_phys(ndim),max_phys(ndim),len_phys(ndim),&
            &         ghostsize(ndim),nm(ndim),nm_ref(ndim),&
            &         h(ndim),stat=info)
        
        min_phys(1:ndim) = 0.0_mk
        max_phys(1:ndim) = 1.0_mk
        len_phys(1:ndim) = max_phys-min_phys
        ghostsize(1:ndim) = 2
        bcdef(1:2*ndim) = ppm_param_bcdef_periodic
        tolexp = -12
        np = 0 
        nullify(field,field_ref)
        nullify(isublist,istart,ndata,istart_ref,ndata_ref)

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

        deallocate(min_phys,max_phys,len_phys,ghostsize,nm)

    end finalize

!----------------------------------------------

!------------- setup --------------------------

    setup

!        call random_seed(size=seedsize)
!        allocate(seed(seedsize))
!        allocate(randnb((1+ndim)*np),stat=info)
!        do i=1,seedsize
!            seed(i)=10+i*i*(rank+1)
!        enddo
!        call random_seed(put=seed)
!        call random_number(randnb)
        

    end setup
!----------------------------------------------
        

!--------------- teardown ---------------------
    teardown
        

    end teardown
!----------------------------------------------


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
