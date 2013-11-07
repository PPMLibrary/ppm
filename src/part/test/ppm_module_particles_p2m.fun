test_suite ppm_module_particles_p2m

use ppm_module_mesh_typedef
use ppm_module_topo_typedef
use ppm_module_field_typedef
use ppm_module_particles_typedef
use ppm_module_mktopo
use ppm_module_io_vtk

#ifdef __MPI
    INCLUDE "mpif.h"
#endif

integer, parameter              :: debug = 0
integer, parameter              :: mk = kind(1.0d0) !kind(1.0e0)
real(mk),parameter              :: pi = ACOS(-1._mk)
integer,parameter               :: ndim=2
integer                         :: decomp,assig,tolexp
integer                         :: info,comm,rank,nproc
real(mk)                        :: tol
integer                         :: topoid=-1
real(mk),dimension(:  ),pointer :: min_phys => NULL()
real(mk),dimension(:  ),pointer :: max_phys => NULL()
real(ppm_kind_double)           :: t2,t1

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
real(mk),dimension(2*ndim)       :: my_patch
real(mk),dimension(ndim)         :: offset

real(mk), dimension(:,:), pointer              :: wp_2r => NULL()


class(ppm_t_main_abstr),POINTER  :: abstr_point => NULL()
TYPE(ppm_v_main_abstr)  :: LFields

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

    test part_to_meshinterp
        type(ppm_t_field) ,TARGET :: VField1,VField2,VField3,VField4
        type(ppm_t_field) ,TARGET :: SField1,SField2,SField3,Vol
        type(ppm_t_particles_d) :: Part1
        real(ppm_kind_double),dimension(ndim) :: pos
        real(ppm_kind_double),dimension(ndim) :: cutoff
        real(ppm_kind_double)                 :: voln
        integer :: np_global

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
        Nm = 31
        np_global = PRODUCT(Nm(1:ndim)-1)
        !Note:
        ! We get an exact p2m interpolation of 2nd order polynomials
        ! (with Mp4) only if the particles are on Cartesian grid with
        ! the same spacing as the mesh).
        call Mesh1%create(topoid,offset,info,Nm=Nm,&
            ghostsize=ighostsize,name='Test_Mesh_1')
            Assert_Equal(info,0)
        !----------------
        ! Add a patch
        !----------------
        if (ndim.eq.2) then
            my_patch(1:4) = (/0.15_mk,0.10_mk,0.99_mk,0.7_mk/)
            my_patch(1:4) = (/0.15_mk,0.15_mk,0.7_mk,0.7_mk/)
        else
            my_patch(1:6) = (/0.15_mk,0.10_mk,0.25_mk,0.89_mk,0.7_mk,0.78_mk/)
        endif
        call Mesh1%def_patch(my_patch,info)
        Assert_Equal(info,0)

        !----------------
        ! Create particles, from a grid + small random displacement
        !----------------
        call Part1%initialize(np_global,info,topoid=topoid,name="Part1")
            Assert_Equal(info,0)

        allocate(wp_2r(ndim,Part1%Npart))
!        call random_number(wp_2r)
        wp_2r(1,:) = -0.1100010000_mk !(wp_2r - 0.5_mk) * Part1%h_avg * 0.0015_mk
        wp_2r(2,:) = 0.1300010000_mk !(wp_2r - 0.5_mk) * Part1%h_avg * 0.0015_mk
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
        call VField1%discretize_on(Part1,info)
        call VField2%create(3,info,name='vecField2')
        call VField2%discretize_on(Part1,info)
        call VField3%create(4,info,name='vecField3')
        call VField3%discretize_on(Part1,info)
        call VField4%create(5,info,name='vecField4')
        call VField4%discretize_on(Part1,info)

        call SField1%create(1,info,name='scaField1')
        call SField1%discretize_on(Part1,info)
        call SField2%create(1,info,name='scaField2')
        call SField2%discretize_on(Part1,info)
        call SField3%create(1,info,name='scaField3')
        call SField3%discretize_on(Part1,info)
        call Vol%create(1,info,name='Part_Volume')
        call Vol%discretize_on(Part1,info)

        !----------------
        ! Initialize the fields with test functions (polynomials of orders 0,1
        ! and 2), to test interpolants of orders up to 3.
        !----------------
        foreach p in particles(Part1) with positions(x) vec_fields(V1=VField1,V2=VField2,V3=VField3,V4=Vfield4) sca_fields(S1=SField1,S2=SField2,S3=SField3,Vol=Vol)
                Vol_p   = 1._mk / np_global
                S1_p    = f_cst(x_p(1:ndim),ndim) * Vol_p
                S2_p    = f_lin(x_p(1:ndim),ndim) * Vol_p
                S3_p    = f_sq(x_p(1:ndim),ndim) * Vol_p

                V1_p(1) = f_cst(x_p(1:ndim),ndim) * Vol_p
                V1_p(2) = f_sq(x_p(1:ndim),ndim) * Vol_p

                V2_p(1) = f_cst(x_p(1:ndim),ndim) * Vol_p
                V2_p(2) = f_lin(x_p(1:ndim),ndim) * Vol_p
                V2_p(3) = f_sq(x_p(1:ndim),ndim) * Vol_p

                V3_p(1) = f_cst(x_p(1:ndim),ndim) * Vol_p
                V3_p(2) = f_lin(x_p(1:ndim),ndim) * Vol_p
                V3_p(3) = f_sq(x_p(1:ndim),ndim) * Vol_p
                V3_p(4) = f_sq(x_p(1:ndim),ndim) * Vol_p

                V4_p(1) = f_cst(x_p(1:ndim),ndim) * Vol_p
                V4_p(2) = f_lin(x_p(1:ndim),ndim) * Vol_p
                V4_p(3) = f_sq(x_p(1:ndim),ndim) * Vol_p
                V4_p(4) = f_sq(x_p(1:ndim),ndim) * Vol_p
                V4_p(5) = f_sq(x_p(1:ndim),ndim) * Vol_p
        end foreach

        !----------------
        ! Get ghost values for all the fields
        !----------------
        call Part1%map_ghosts(info)
        Assert_Equal(info,0)


        abstr_point => SField2
        CALL LFields%push(abstr_point,info)
        Assert_Equal(info,0)
        CALL ppm_vtk_particles("output",Part1,info,Fields=LFields)
        Assert_Equal(info,0)
        abstr_point => SField1
        CALL LFields%push(abstr_point,info)
        Assert_Equal(info,0)
        CALL ppm_vtk_particles("output",Part1,info,Fields=LFields)
        Assert_Equal(info,0)
        abstr_point => VField1
        CALL LFields%push(abstr_point,info)
        CALL ppm_vtk_particles("output",Part1,info,Fields=LFields)
        Assert_Equal(info,0)
        abstr_point => VField2
        CALL LFields%push(abstr_point,info)
        CALL ppm_vtk_particles("output",Part1,info,Fields=LFields)
        Assert_Equal(info,0)
        abstr_point => VField4
        CALL LFields%push(abstr_point,info)
        abstr_point => VField3
        CALL LFields%push(abstr_point,info)
        abstr_point => SField3
        CALL LFields%push(abstr_point,info)
        CALL ppm_vtk_particles("output",Part1,info,Fields=LFields)
        Assert_Equal(info,0)
        CALL ppm_vtk_particles("output",Part1,info)
        Assert_Equal(info,0)
        !----------------
        ! Perform the p2m interpolation
        !----------------

        CALL ppm_vtk_particles("output_before",Part1,info)
        Assert_Equal(info,0)

        call Part1%interp_to_mesh(Mesh1,VField1,ppm_param_rmsh_kernel_mp4,info)
            Assert_Equal(info,0)
        call Part1%interp_to_mesh(Mesh1,VField2,ppm_param_rmsh_kernel_mp4,info)
            Assert_Equal(info,0)
        call Part1%interp_to_mesh(Mesh1,VField3,ppm_param_rmsh_kernel_mp4,info)
            Assert_Equal(info,0)
        call Part1%interp_to_mesh(Mesh1,VField4,ppm_param_rmsh_kernel_mp4,info)
            Assert_Equal(info,0)

        call Part1%interp_to_mesh(Mesh1,SField1,ppm_param_rmsh_kernel_mp4,info)
            Assert_Equal(info,0)
        call Part1%interp_to_mesh(Mesh1,SField2,ppm_param_rmsh_kernel_mp4,info)
            Assert_Equal(info,0)
        call Part1%interp_to_mesh(Mesh1,SField3,ppm_param_rmsh_kernel_mp4,info)
            Assert_Equal(info,0)

        CALL ppm_vtk_particles("output_after",Part1,info)
        Assert_Equal(info,0)


        tol = 1e-2

        !----------------
        ! Define a cutoff distance from the sides of the patch. Particles that
        ! are in the patch but too close to the sides will not receive anything
        ! from the patch during interpolation (except if there are some periodic
        ! boundaries, of course...)
        !----------------
        cutoff = REAL(Mesh1%ghostsize(1:ndim),ppm_kind_double)* Mesh1%h(1:ndim)

        voln = 1._mk / PRODUCT(Mesh1%Nm(1:ndim)-1)

        !----------------
        !Loop through all the mesh nodes and check that the values of the field
        !have been interpolated EXACTLY,
        !as should be the case for constant, linear and quadratic functions
        !with a 3rd order interpolating kernel (like Mp4)
        !----------------
        IF (ppm_debug.GE.1) THEN
        IF (ndim.EQ.2) THEN
        foreach n in equi_mesh(Mesh1) with sca_fields(SField1,SField2,SField3) indices(i,j)
            for real
                pos(1:ndim) = sbpitr%get_pos(i,j)
                write(100+rank,'(9(E12.5,1X))') pos(1),pos(2),pos(ndim),&
                    f_cst(pos(1:ndim),ndim), &
                    f_lin(pos(1:ndim),ndim),&
                    f_sq(pos(1:ndim),ndim),&
                    SField1_n/voln,&
                    SField2_n/voln,&
                    SField3_n/voln
        end foreach
        ELSE
        foreach n in equi_mesh(Mesh1) with sca_fields(SField1,SField2,SField3) indices(i,j,k)
            for real
                pos(1:ndim) = sbpitr%get_pos(i,j,k)
                write(100+rank,'(9(E12.5,1X))') pos(1),pos(2),pos(ndim),&
                    f_cst(pos(1:ndim),ndim), &
                    f_lin(pos(1:ndim),ndim),&
                    f_sq(pos(1:ndim),ndim),&
                    SField1_n/voln,&
                    SField2_n/voln,&
                    SField3_n/voln
        end foreach
        ENDIF

        foreach p in particles(Part1) with positions(x) sca_fields(S1=SField1,S2=SField2,S3=SField3,V=Vol)
                write(200+rank,'(6(E12.5,1X))') x_p(1),x_p(2),x_p(ndim),&
                    S1_p/V_p, S2_p/V_p, S3_p/V_p
        end foreach
        ENDIF
        !----------------
        ! Check what we've got
        !----------------
        IF (ndim.EQ.2) THEN
        foreach n in equi_mesh(Mesh1) with sca_fields(SField1,SField2,SField3) vec_fields(VField1,VField2,VField3,VField4) indices(i,j)
            for real
                pos(1:ndim) = sbpitr%get_pos(i,j)
                Assert_Equal_Within(VField1_n(1)/voln,f_cst(pos(1:ndim),ndim),tol)
                Assert_Equal_Within(VField1_n(2)/voln,f_sq(pos(1:ndim),ndim),tol)

                Assert_Equal_Within(VField2_n(1)/voln,f_cst(pos(1:ndim),ndim),tol)
                Assert_Equal_Within(VField2_n(2)/voln,f_lin(pos(1:ndim),ndim),tol)
                Assert_Equal_Within(VField2_n(3)/voln,f_sq(pos(1:ndim),ndim),tol)

                Assert_Equal_Within(VField3_n(1)/voln,f_cst(pos(1:ndim),ndim),tol)
                Assert_Equal_Within(VField3_n(2)/voln,f_lin(pos(1:ndim),ndim),tol)
                Assert_Equal_Within(VField3_n(3)/voln,f_sq(pos(1:ndim),ndim),tol)
                Assert_Equal_Within(VField3_n(4)/voln,f_sq(pos(1:ndim),ndim),tol)

                Assert_Equal_Within(VField4_n(1)/voln,f_cst(pos(1:ndim),ndim),tol)
                Assert_Equal_Within(VField4_n(2)/voln,f_lin(pos(1:ndim),ndim),tol)
                Assert_Equal_Within(VField4_n(3)/voln,f_sq(pos(1:ndim),ndim),tol)
                Assert_Equal_Within(VField4_n(4)/voln,f_sq(pos(1:ndim),ndim),tol)
                Assert_Equal_Within(VField4_n(5)/voln,f_sq(pos(1:ndim),ndim),tol)

                Assert_Equal_Within(SField1_n/voln,f_cst(pos(1:ndim),ndim),tol)
                Assert_Equal_Within(SField2_n/voln,f_lin(pos(1:ndim),ndim),tol)
                Assert_Equal_Within(SField3_n/voln,f_sq(pos(1:ndim),ndim),tol)
        end foreach

        ELSE

        !foreach n in equi_mesh(Mesh1) with sca_fields(SField1,SField2,SField3) vec_fields(VField1,VField2,VField3,VField4) indices(i,j,k)
        foreach n in equi_mesh(Mesh1) with sca_fields(SField1,SField2,SField3) vec_fields(VField1,VField2,VField3,VField4) indices(i,j,k)
            for real
                pos(1:ndim) = sbpitr%get_pos(i,j,k)
                Assert_Equal_Within(VField1_n(1)/voln,f_cst(pos(1:ndim),ndim),tol)
                Assert_Equal_Within(VField1_n(2)/voln,f_sq(pos(1:ndim),ndim),tol)

                Assert_Equal_Within(VField2_n(1)/voln,f_cst(pos(1:ndim),ndim),tol)
                Assert_Equal_Within(VField2_n(2)/voln,f_lin(pos(1:ndim),ndim),tol)
                Assert_Equal_Within(VField2_n(3)/voln,f_sq(pos(1:ndim),ndim),tol)

                Assert_Equal_Within(VField3_n(1)/voln,f_cst(pos(1:ndim),ndim),tol)
                Assert_Equal_Within(VField3_n(2)/voln,f_lin(pos(1:ndim),ndim),tol)
                Assert_Equal_Within(VField3_n(3)/voln,f_sq(pos(1:ndim),ndim),tol)
                Assert_Equal_Within(VField3_n(4)/voln,f_sq(pos(1:ndim),ndim),tol)

                Assert_Equal_Within(VField4_n(1)/voln,f_cst(pos(1:ndim),ndim),tol)
                Assert_Equal_Within(VField4_n(2)/voln,f_lin(pos(1:ndim),ndim),tol)
                Assert_Equal_Within(VField4_n(3)/voln,f_sq(pos(1:ndim),ndim),tol)
                Assert_Equal_Within(VField4_n(4)/voln,f_sq(pos(1:ndim),ndim),tol)
                Assert_Equal_Within(VField4_n(5)/voln,f_sq(pos(1:ndim),ndim),tol)

                Assert_Equal_Within(SField1_n/voln,f_cst(pos(1:ndim),ndim),tol)
                Assert_Equal_Within(SField2_n/voln,f_lin(pos(1:ndim),ndim),tol)
                Assert_Equal_Within(SField3_n/voln,f_sq(pos(1:ndim),ndim),tol)
        end foreach

        ENDIF


        CALL MPI_BARRIER(comm,info)


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
        !check that we are leaving the test without error
        Assert_Equal(info,0)
    end test

!-------------------------------------------------------------
! test function
!-------------------------------------------------------------
pure function f_cst(pos,ndim) RESULT(res)
    real(mk)                              :: res
    integer                 ,  intent(in) :: ndim
    real(mk), dimension(ndim), intent(in) :: pos

    res =  1._mk !42.17_mk
end function

pure function f_lin(pos,ndim) RESULT(res)
    real(mk)                              :: res
    integer                 ,  intent(in) :: ndim
    real(mk), dimension(ndim), intent(in) :: pos

    res =  pos(1) + 10._mk*pos(2) + 100._mk*pos(ndim)
end function

pure function f_sq(pos,ndim) RESULT(res)
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
