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
real(mk),dimension(2*ndim)       :: my_patch
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
        type(ppm_t_field) :: Field1,Field2,Field3,Field4
        type(ppm_t_particles_d) :: Part1

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


        Nm(1) = 35
        Nm(ndim) = 24
        Nm(2) = 65
        call Mesh1%create(topoid,offset,info,Nm=Nm,&
            ghostsize=ighostsize,name='Test_Mesh_1')
            Assert_Equal(info,0)

        call Part1%initialize(np_global,info,topoid=topoid,name="Part1")
            Assert_Equal(info,0)

        call Part1%set_cutoff(3._mk * Part1%h_avg,info)
        Assert_Equal(info,0)

        allocate(wp_2r(ndim,Part1%Npart))
        call random_number(wp_2r)
        wp_2r = (wp_2r - 0.5_mk) * Part1%h_avg * 0.15_mk
        call Part1%move(wp_2r,info)
        Assert_Equal(info,0)
        deallocate(wp_2r)

        call Part1%apply_bc(info)
        Assert_Equal(info,0)

        call Part1%map(info,global=.true.,topoid=topoid)
        Assert_Equal(info,0)

        if (ndim.eq.2) then
            my_patch(1:4) = (/0.15_mk,0.10_mk,0.99_mk,0.7_mk/)
        else
            my_patch(1:6) = (/0.15_mk,0.10_mk,0.5_mk,0.89_mk,0.7_mk,0.78_mk/)
        endif
        call Mesh1%def_patch(my_patch,info) 
        Assert_Equal(info,0)

        call Field1%create(3,info,name='vecField') 
        call Field1%discretize_on(Mesh1,info)

        call Field2%create(1,info,name='scaField1') 
        call Field2%discretize_on(Mesh1,info)

        call Field3%create(1,info,name='scaField2') 
        call Field3%discretize_on(Mesh1,info)

        call Field4%create(1,info,name='scaField3') 
        call Field4%discretize_on(Mesh1,info)

        IF (ndim.EQ.2) THEN
        foreach n in equi_mesh(Mesh1) with sca_fields(Field2,Field3,Field4) vec_fields(Field1) indices(i,j)
            for all
                pos(1:ndim) = sbpitr%get_pos(i,j)
                Field1_n(1) = test_constant(pos(1:ndim),ndim)
                Field1_n(2) = test_linear(pos(1:ndim),ndim)
                Field1_n(3) = test_quadratic(pos(1:ndim),ndim)
                Field2_n    = test_constant(pos(1:ndim),ndim)
                Field3_n    = test_linear(pos(1:ndim),ndim)
                Field4_n    = test_quadratic(pos(1:ndim),ndim)
        end foreach

        ELSE

        foreach n in equi_mesh(Mesh1) with sca_fields(Field2,Field3,Field4) vec_fields(Field1) indices(i,j,k)
            for all
                pos(1:ndim) = sbpitr%get_pos(i,j,k)
                Field1_n(1) = test_constant(pos(1:ndim),ndim)
                Field1_n(2) = test_linear(pos(1:ndim),ndim)
                Field1_n(3) = test_quadratic(pos(1:ndim),ndim)
                Field2_n    = test_constant(pos(1:ndim),ndim)
                Field3_n    = test_linear(pos(1:ndim),ndim)
                Field4_n    = test_quadratic(pos(1:ndim),ndim)
        end foreach

        ENDIF

        call Mesh1%interp_to_part(Part1,Field1,ppm_param_rmsh_kernel_mp4,info)
            Assert_Equal(info,0)
        call Mesh1%interp_to_part(Part1,Field2,ppm_param_rmsh_kernel_mp4,info)
            Assert_Equal(info,0)
        call Mesh1%interp_to_part(Part1,Field3,ppm_param_rmsh_kernel_mp4,info)
            Assert_Equal(info,0)
        call Mesh1%interp_to_part(Part1,Field4,ppm_param_rmsh_kernel_mp4,info)
            Assert_Equal(info,0)

        tol = 1e-12

        cutoff = REAL(Mesh1%ghostsize(1:ndim),ppm_kind_double)* Mesh1%h(1:ndim)

        !Loop through all particles and check that the values of the field
        !have been interpolated EXACTLY,
        !as should be the case for constant, linear and quadratic functions
        !with a 3rd order interpolating kernel (like Mp4)
        foreach p in particles(Part1) with positions(x) vec_fields(u=Field1) sca_fields(sca1=Field2,sca2=Field3,sca3=Field4)
            if (is_well_within(x_p(1:ndim),my_patch(1:2*ndim),cutoff,ndim)) then
! Lots of debugging output in case something goes wrong....
!                IF (abs(u_p(3)-test_quadratic(x_p(1:ndim),ndim)).GT.tol) then
!                    !find to which subpatch this particle belongs
!                    patch => Mesh1%subpatch%begin()
!                    ploop: do while(associated(patch))
!                     IF( (      x_p(1).GE.patch%start_red(1) .AND. &
!                     &          x_p(2).GE.patch%start_red(2) .AND. &
!                     &          x_p(3).GE.patch%start_red(3) .AND. &
!                     &          x_p(1).LE.patch%end_red(1) .AND. &
!                     &          x_p(3).LE.patch%end_red(3) .AND. &
!                     &          x_p(2).LE.patch%end_red(2) ) ) THEN
!                                          IF(   (x_p(1).LT.patch%end(1) .OR.  &
!                     &                           (patch%bc(2).GE.0   .AND.    &
!                     &                           patch%bc(2).NE. ppm_param_bcdef_periodic)).AND.&
!                     &                          (x_p(2).LT.patch%end(2) .OR.  &
!                     &                           (patch%bc(4).GE.0   .AND.    &
!                     &                           patch%bc(4).NE. ppm_param_bcdef_periodic)).AND.&
!                     &                          (x_p(3).LT.patch%end(3) .OR.  &
!                     &                           (patch%bc(6).GE.0   .AND.    &
!                     &                           patch%bc(6).NE. ppm_param_bcdef_periodic))) THEN
!                                          EXIT ploop
!                                          ENDIF
!                    ENDIF
!                    patch => Mesh1%subpatch%next()
!                    enddo ploop
!                    check_associated(patch,"particle does not belong to any subpatch")
!                    stdout_f('(A,3(F17.13,1X))',"part pos ",'x_p(1:ndim)')
!                    stdout_f('(A,6(F7.3,1X))',"   is inside patch",&
!                        'my_patch(1:2*ndim)')
!                    stdout_f('(A,6(F7.3,1X))',"patch%start     = ",'patch%start')
!                    stdout_f('(A,6(F7.3,1X))',"patch%end       = ",'patch%end')
!                    stdout_f('(A,6(I7  ,1X))',"patch%istart    = ",'patch%istart')
!                    stdout_f('(A,6(I7  ,1X))',"patch%iend      = ",'patch%iend')
!                    stdout_f('(A,6(F7.3,1X))',"patch%istart    = ",&
!                        '(patch%istart-1) * Mesh1%h')
!                    stdout_f('(A,6(F7.3,1X))',"patch%iend      = ",&
!                        '(patch%iend-1) * Mesh1%h')
!                    stdout_f('(A,3(I7,1X))',  "patch%nnodes    = ",'patch%nnodes')
!                    stdout_f('(A,6(F7.3,1X))',"patch%start_red = ",'patch%start_red')
!                    stdout_f('(A,6(F7.3,1X))',"patch%end_red   = ",'patch%end_red')
!                    stdout("patch%ghostsize",'patch%ghostsize')
!                    x0 = x_p(1:ndim)/Mesh1%h-patch%istart+1
!                    stdout_f('(A,3(I2,1X))',"ipj0= ",'FLOOR(x0)')
!                    stdout_f('(A,3(I2,1X))',"ipj1= ",'FLOOR(x0)+1')
!                    stdout_f('(A,3(I2,1X))',"ipj2= ",'FLOOR(x0)+2')
!                    stdout_f('(A,3(I2,1X))',"ipj3= ",'FLOOR(x0)+3')
!                    stdout_f('(A,3(F7.3,1X))',"xpj= ",'x0-REAL(FLOOR(x0))')
!                    stdout_f('(A,3(E13.4,1X))',"Mesh1%h",'Mesh1%h')
!                ENDIF

                Assert_Equal_Within(sca1_p,test_constant(x_p(1:ndim),ndim),tol)
                Assert_Equal_Within(sca2_p,test_linear(x_p(1:ndim),ndim),tol)
                Assert_Equal_Within(sca3_p,test_quadratic(x_p(1:ndim),ndim),tol)
                Assert_Equal_Within(u_p(1),test_constant(x_p(1:ndim),ndim),tol)
                Assert_Equal_Within(u_p(2),test_linear(x_p(1:ndim),ndim),tol)
                Assert_Equal_Within(u_p(3),test_quadratic(x_p(1:ndim),ndim),tol)

            endif
        end foreach

        CALL MPI_BARRIER(comm,info)


        call Mesh1%destroy(info)
            Assert_Equal(info,0)
        call Field1%destroy(info)
            Assert_Equal(info,0)
        call Field2%destroy(info)
            Assert_Equal(info,0)
        call Part1%destroy(info)
            Assert_Equal(info,0)

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
