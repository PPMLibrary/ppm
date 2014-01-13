test_suite ppm_module_mesh_m2p

use ppm_module_mesh_typedef
use ppm_module_topo_typedef
use ppm_module_field_typedef
use ppm_module_particles_typedef
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
integer                         :: seedsize
integer,  dimension(:),allocatable :: seed

integer                         :: i,j,k
integer                         :: nsublist
integer, dimension(:  ),pointer :: isublist => NULL()
integer, dimension(2*ndim)      :: bcdef
real(mk),dimension(:  ),pointer :: cost => NULL()
integer, dimension(:  ),pointer :: nm => NULL()
real(mk),dimension(:  ),pointer :: h => NULL()
type(ppm_t_topo),       pointer :: topo => NULL()

type(ppm_t_equi_mesh),TARGET     :: Mesh1,Mesh2
integer                          :: ipatch,isub,jsub
class(ppm_t_subpatch_),POINTER   :: p => NULL()
class(ppm_t_subpatch_),POINTER   :: patch => NULL()

integer                          :: mypatchid
real(mk),dimension(6)       :: my_patch
real(mk),dimension(ndim)         :: offset

real(mk), dimension(:,:), pointer              :: wp_2r => NULL()

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

        call random_seed(size=seedsize)
        allocate(seed(seedsize))
        do i=1,seedsize
            seed(i)=9+i*(rank+1)
        enddo
        call random_seed(put=seed)

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

    test mesh_to_part_interp
        type(ppm_t_field) :: VField1,VField2,VField3,VField4
        type(ppm_t_field) :: SField1,SField2,SField3
        type(ppm_t_particles_d),TARGET :: Part1

        real(ppm_kind_double),dimension(ndim) :: x0

        real(ppm_kind_double),dimension(ndim) :: pos
        real(ppm_kind_double),dimension(ndim) :: cutoff
        integer :: np_global = 500000

        start_subroutine("test_interp")

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        assig  = ppm_param_assign_internal
        topoid = 0
        sca_ghostsize = 0.07_mk
        call ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,    &
            &               bcdef,sca_ghostsize,cost,info)
        Assert_Equal(info,0)


        !----------------
        ! Create Mesh
        !----------------
        Nm(1) = 35
        Nm(ndim) = 24
        Nm(2) = 65
        call Mesh1%create(topoid,offset,info,Nm=Nm,&
            ghostsize=ighostsize,name='Test_Mesh_1')
            Assert_Equal(info,0)
        !----------------
        ! Add a patch
        !----------------
        if (ndim.eq.2) then
            my_patch(1:4) = (/0.15_mk,0.10_mk,0.99_mk,0.7_mk/)
        else
            my_patch(1:6) = (/0.15_mk,0.10_mk,0.5_mk,0.89_mk,0.7_mk,0.78_mk/)
        endif
        call Mesh1%def_patch(my_patch,info)
        Assert_Equal(info,0)

        !----------------
        ! Create particles, from a grid + small random displacement
        !----------------
        call Part1%initialize(np_global,info,topoid=topoid,name="Part1")
            Assert_Equal(info,0)

        allocate(wp_2r(ndim,Part1%Npart))
        call random_number(wp_2r)
        wp_2r = (wp_2r - 0.5_mk) * Part1%h_avg * 1.15_mk
        call Part1%move(wp_2r,info)
        Assert_Equal(info,0)
        deallocate(wp_2r)

        !----------------
        ! Put particles back into the domain and global map them
        !----------------
        call Part1%apply_bc(info)
        Assert_Equal(info,0)

        call Part1%map(info,global=.true.,topoid=topoid)
        Assert_Equal(info,0)


        !----------------
        ! Define some fields. Vector and scalar fields, with different
        ! dimensions because the interpolation routines are hard-coded for some
        ! and we want to test them all!
        !----------------
        call VField1%create(2,info,name='vecField1')
        call VField1%discretize_on(Mesh1,info)
        call VField2%create(3,info,name='vecField2')
        call VField2%discretize_on(Mesh1,info)
        call VField3%create(4,info,name='vecField3')
        call VField3%discretize_on(Mesh1,info)
        call VField4%create(5,info,name='vecField4')
        call VField4%discretize_on(Mesh1,info)

        call SField1%create(1,info,name='scaField1')
        call SField1%discretize_on(Mesh1,info)
        call SField2%create(1,info,name='scaField2')
        call SField2%discretize_on(Mesh1,info)
        call SField3%create(1,info,name='scaField3')
        call SField3%discretize_on(Mesh1,info)

        !----------------
        ! Initialize the fields with test functions (polynomials of orders 0,1
        ! and 2), to test interpolants of orders up to 3.
        !----------------
        IF (ndim.EQ.2) THEN
        foreach n in equi_mesh(Mesh1) with sca_fields(SField1,SField2,SField3) vec_fields(VField1,VField2,VField3,VField4) indices(i,j)
            for all
                pos(1:ndim) = sbpitr%get_pos(i,j)
                VField1_n(1) = test_constant(pos(1:ndim),ndim)
                VField1_n(2) = test_quadratic(pos(1:ndim),ndim)

                VField2_n(1) = test_constant(pos(1:ndim),ndim)
                VField2_n(2) = test_linear(pos(1:ndim),ndim)
                VField2_n(3) = test_quadratic(pos(1:ndim),ndim)

                VField3_n(1) = test_constant(pos(1:ndim),ndim)
                VField3_n(2) = test_linear(pos(1:ndim),ndim)
                VField3_n(3:4) = test_quadratic(pos(1:ndim),ndim)

                VField4_n(1) = test_constant(pos(1:ndim),ndim)
                VField4_n(2) = test_linear(pos(1:ndim),ndim)
                VField4_n(3:5) = test_quadratic(pos(1:ndim),ndim)

                SField1_n    = test_constant(pos(1:ndim),ndim)
                SField2_n    = test_linear(pos(1:ndim),ndim)
                SField3_n    = test_quadratic(pos(1:ndim),ndim)
        end foreach

        ELSE

        foreach n in equi_mesh(Mesh1) with sca_fields(SField1,SField2,SField3) vec_fields(VField1,VField2,VField3,VField4) indices(i,j,k)
            for all
                pos(1:ndim) = sbpitr%get_pos(i,j,k)
                VField1_n(1) = test_constant(pos(1:ndim),ndim)
                VField1_n(2) = test_quadratic(pos(1:ndim),ndim)

                VField2_n(1) = test_constant(pos(1:ndim),ndim)
                VField2_n(2) = test_linear(pos(1:ndim),ndim)
                VField2_n(3) = test_quadratic(pos(1:ndim),ndim)

                VField3_n(1) = test_constant(pos(1:ndim),ndim)
                VField3_n(2) = test_linear(pos(1:ndim),ndim)
                VField3_n(3:4) = test_quadratic(pos(1:ndim),ndim)

                VField4_n(1) = test_constant(pos(1:ndim),ndim)
                VField4_n(2) = test_linear(pos(1:ndim),ndim)
                VField4_n(3:5) = test_quadratic(pos(1:ndim),ndim)

                SField1_n    = test_constant(pos(1:ndim),ndim)
                SField2_n    = test_linear(pos(1:ndim),ndim)
                SField3_n    = test_quadratic(pos(1:ndim),ndim)
        end foreach

        ENDIF

        !----------------
        ! Perform the m2p interpolation
        !----------------
        call Mesh1%interp_to_part(Part1,VField1,ppm_param_rmsh_kernel_mp4,info)
            Assert_Equal(info,0)
        call Mesh1%interp_to_part(Part1,VField2,ppm_param_rmsh_kernel_mp4,info)
            Assert_Equal(info,0)
        call Mesh1%interp_to_part(Part1,VField3,ppm_param_rmsh_kernel_mp4,info)
            Assert_Equal(info,0)
        call Mesh1%interp_to_part(Part1,VField4,ppm_param_rmsh_kernel_mp4,info)
            Assert_Equal(info,0)
        call Mesh1%interp_to_part(Part1,SField1,ppm_param_rmsh_kernel_mp4,info)
            Assert_Equal(info,0)
        call Mesh1%interp_to_part(Part1,SField2,ppm_param_rmsh_kernel_mp4,info)
            Assert_Equal(info,0)
        call Mesh1%interp_to_part(Part1,SField3,ppm_param_rmsh_kernel_mp4,info)
            Assert_Equal(info,0)

        tol = 1e-12

        !----------------
        ! Define a cutoff distance from the sides of the patch. Particles that
        ! are in the patch but too close to the sides will not receive anything
        ! from the patch during interpolation (except if there are some periodic
        ! boundaries, of course...)
        !----------------
        cutoff = REAL(Mesh1%ghostsize(1:ndim),ppm_kind_double)* Mesh1%h(1:ndim)

        !----------------
        !Loop through all particles and check that the values of the field
        !have been interpolated EXACTLY,
        !as should be the case for constant, linear and quadratic functions
        !with a 3rd order interpolating kernel (like Mp4)
        !----------------
        foreach p in particles(Part1) with positions(x) vec_fields(V1=VField1,V2=VField2,V3=VField3,V4=Vfield4) sca_fields(S1=SField1,S2=SField2,S3=SField3)
            if (is_well_within(x_p(1:ndim),my_patch(1:2*ndim),cutoff,ndim)) then

                Assert_Equal_Within(S1_p,test_constant(x_p(1:ndim),ndim),tol)
                Assert_Equal_Within(S2_p,test_linear(x_p(1:ndim),ndim),tol)
                Assert_Equal_Within(S3_p,test_quadratic(x_p(1:ndim),ndim),tol)

                Assert_Equal_Within(V1_p(1),test_constant(x_p(1:ndim),ndim),tol)
                Assert_Equal_Within(V1_p(2),test_quadratic(x_p(1:ndim),ndim),tol)

                Assert_Equal_Within(V2_p(1),test_constant(x_p(1:ndim),ndim),tol)
                Assert_Equal_Within(V2_p(2),test_linear(x_p(1:ndim),ndim),tol)
                Assert_Equal_Within(V2_p(3),test_quadratic(x_p(1:ndim),ndim),tol)

                Assert_Equal_Within(V3_p(1),test_constant(x_p(1:ndim),ndim),tol)
                Assert_Equal_Within(V3_p(2),test_linear(x_p(1:ndim),ndim),tol)
                Assert_Equal_Within(V3_p(3),test_quadratic(x_p(1:ndim),ndim),tol)
                Assert_Equal_Within(V3_p(4),test_quadratic(x_p(1:ndim),ndim),tol)

                Assert_Equal_Within(V4_p(1),test_constant(x_p(1:ndim),ndim),tol)
                Assert_Equal_Within(V4_p(2),test_linear(x_p(1:ndim),ndim),tol)
                Assert_Equal_Within(V4_p(3),test_quadratic(x_p(1:ndim),ndim),tol)
                Assert_Equal_Within(V4_p(4),test_quadratic(x_p(1:ndim),ndim),tol)
                Assert_Equal_Within(V4_p(5),test_quadratic(x_p(1:ndim),ndim),tol)

            endif
        end foreach
#ifdef __MPI
        CALL MPI_BARRIER(comm,info)
#endif

        call Mesh1%destroy(info)
        call SField1%destroy(info)
        call SField2%destroy(info)
        call SField3%destroy(info)
        call VField1%destroy(info)
        call VField2%destroy(info)
        call VField3%destroy(info)
        call VField4%destroy(info)
        call Part1%destroy(info)

        end_subroutine()
    end test

!-------------------------------------------------------------
! test function
!-------------------------------------------------------------
pure function test_constant(pos,ndim) RESULT(res)
    real(mk)                              :: res
    integer                 ,  intent(in) :: ndim
    real(mk), dimension(ndim), intent(in) :: pos

    res =  42.17_mk
end function

pure function test_linear(pos,ndim) RESULT(res)
    real(mk)                              :: res
    integer                 ,  intent(in) :: ndim
    real(mk), dimension(ndim), intent(in) :: pos

    res =  pos(1) + 10._mk*pos(2) + 100._mk*pos(ndim)
end function

pure function test_quadratic(pos,ndim) RESULT(res)
    real(mk)                              :: res
    integer                 ,  intent(in) :: ndim
    real(mk), dimension(ndim), intent(in) :: pos

    res =  pos(1)**2 + 10._mk*pos(2)**2 + 100._mk*pos(ndim)**2
end function

!!! check whether a particle is within a patch and more than a cutoff
!!! distance away from its boundaries.
pure function is_well_within(pos,patch,cutoff,ndim) RESULT(res)
    logical                               :: res
    real(mk), dimension(ndim), intent(in) :: pos
    real(mk), dimension(2*ndim),intent(in):: patch
    real(mk), dimension(ndim), intent(in) :: cutoff
    integer                 ,  intent(in) :: ndim

    res = ALL(pos(1:ndim).GE.(patch(1:ndim)+cutoff(1:ndim)))
    res = res .AND. ALL(pos(1:ndim).LE.(patch(ndim+1:2*ndim)-cutoff(1:ndim)))

end function




end test_suite
