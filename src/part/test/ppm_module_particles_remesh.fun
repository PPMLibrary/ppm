test_suite ppm_module_particles_remesh

USE ppm_module_substart
USE ppm_module_substop
USE ppm_module_mesh_typedef
USE ppm_module_topo_typedef
USE ppm_module_field_typedef
USE ppm_module_particles_typedef
USE ppm_module_mktopo
USE ppm_module_io_vtk

INTEGER, PARAMETER              :: debug = 0
INTEGER, PARAMETER              :: MK = KIND(1.0d0) !KIND(1.0e0)
REAL(MK),PARAMETER              :: pi = ACOS(-1.0_MK)
INTEGER,PARAMETER               :: ndim=2
INTEGER                         :: decomp,assig,tolexp
INTEGER                         :: info,comm,rank,nproc
REAL(MK)                        :: tol
INTEGER                         :: topoid=-1
REAL(MK),DIMENSION(:  ),POINTER :: min_phys => NULL()
REAL(MK),DIMENSION(:  ),POINTER :: max_phys => NULL()
REAL(ppm_kind_double)           :: t2,t1

INTEGER, DIMENSION(:  ),POINTER :: ighostsize => NULL()
REAL(MK)                        :: sca_ghostsize
INTEGER                         :: seedsize
INTEGER,  DIMENSION(:),ALLOCATABLE :: seed

INTEGER                         :: i,j,k
INTEGER                         :: nsublist
INTEGER, DIMENSION(:  ),POINTER :: isublist => NULL()
INTEGER, DIMENSION(2*ndim)      :: bcdef
REAL(MK),DIMENSION(:  ),POINTER :: cost => NULL()
INTEGER, DIMENSION(:  ),POINTER :: nm => NULL()
REAL(MK),DIMENSION(:  ),POINTER :: h => NULL()
TYPE(ppm_t_topo),       POINTER :: topo => NULL()

TYPE(ppm_t_equi_mesh),TARGET     :: Mesh1,Mesh2
INTEGER                          :: ipatch,isub,jsub
CLASS(ppm_t_subpatch_),POINTER   :: p => NULL()
CLASS(ppm_t_subpatch_),POINTER   :: patch => NULL()

INTEGER                          :: mypatchid
REAL(MK),DIMENSION(6)       :: my_patch
REAL(MK),DIMENSION(ndim)         :: offset

REAL(MK), DIMENSION(:,:), POINTER              :: wp_2r => NULL()


CLASS(ppm_t_main_abstr),POINTER  :: abstr_point => NULL()
TYPE(ppm_v_main_abstr)  :: LFields

TYPE(ppm_t_field) ,TARGET :: VField1,VField2,VField3,VField4
TYPE(ppm_t_field) ,TARGET :: SField1,SField2,SField3,Vol
TYPE(ppm_t_particles_d),TARGET :: Part1
REAL(ppm_kind_double),DIMENSION(ndim) :: pos
REAL(ppm_kind_double),DIMENSION(ndim) :: cutoff
REAL(ppm_kind_double)                 :: voln
INTEGER :: np_global

!---------------- init -----------------------

    init

        USE ppm_module_topo_typedef
        USE ppm_module_init

        ALLOCATE(min_phys(ndim),max_phys(ndim),&
            &         ighostsize(ndim),nm(ndim),h(ndim))

        min_phys(1:ndim) = 0.0_MK
        max_phys(1:ndim) = 1.0_MK
        ighostsize(1:ndim) = 2
        bcdef(1:2*ndim) = ppm_param_bcdef_periodic
        tolexp = -12

#ifdef __MPI
        comm = MPI_COMM_WORLD
        CALL MPI_Comm_rank(comm,rank,info)
        CALL MPI_Comm_size(comm,nproc,info)
#else
        rank = 0
        nproc = 1
#endif
        CALL ppm_init(ndim,mk,tolexp,0,debug,info,99)

        CALL RANDOM_SEED(SIZE=seedsize)
        ALLOCATE(seed(seedsize))
        DO i=1,seedsize
            seed(i)=9+i*(rank+1)
        ENDDO
        CALL RANDOM_SEED(PUT=seed)

    end init

!----------------------------------------------

!---------------- finalzie --------------------


    finalize
        USE ppm_module_finalize

        CALL ppm_finalize(info)

        DEALLOCATE(min_phys,max_phys,ighostsize,h,nm)

    end finalize

!----------------------------------------------

!------------- setup --------------------------

    setup
        start_subroutine("setup")

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        assig  = ppm_param_assign_internal
        topoid = 0
        sca_ghostsize = 0.07_MK
        CALL ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,    &
            &               bcdef,sca_ghostsize,cost,info)

        !----------------
        ! Create Mesh
        !----------------
        Nm = 61
        np_global = PRODUCT(Nm(1:ndim)-1)
        !Note:
        ! We get an exact p2m interpolation of 2nd order polynomials
        ! (with Mp4) only if the particles are on Cartesian grid with
        ! the same spacing as the mesh).
        CALL Mesh1%create(topoid,offset,info,Nm=Nm,&
            ghostsize=ighostsize,name='Test_Mesh_1')

        !----------------
        ! Add a patch
        !----------------

!not for now... (first try with a mesh that covers the whole domain)
!               (which is the default, when no patches are defined)

        IF (ndim.EQ.2) THEN
            my_patch(1:4) = (/0.15_MK,0.10_MK,0.99_MK,0.7_MK/)
            my_patch(1:4) = (/0.15_MK,0.15_MK,0.7_MK,0.7_MK/)
            !works (at least on one proc)
            my_patch(1:4) = (/-10000.0_MK,-1000.0_MK,1000.0_MK,1000.0_MK/)
            !does not work (problem with ghosts)
        else
            my_patch(1:6) = (/0.15_MK,0.10_MK,0.25_MK,0.89_MK,0.7_MK,0.78_MK/)
        ENDIF
        CALL Mesh1%def_patch(my_patch,info)

        !----------------
        ! Create particles, from a grid + small random displacement
        !----------------
        CALL Part1%initialize(np_global,info,topoid=topoid,name="Part1")

        ALLOCATE(wp_2r(ndim,Part1%Npart))
!        CALL RANDOM_NUMBER(wp_2r)
!        wp_2r = (wp_2r-0.5_MK)*Part1%h_avg * 0.3_MK
        wp_2r = Part1%h_avg*0.04997_MK

        CALL Part1%move(wp_2r,info)
        DEALLOCATE(wp_2r)

        !----------------
        ! Put particles back into the domain and global map them
        !----------------
        CALL Part1%apply_bc(info)

        CALL Part1%map(info,global=.TRUE.,topoid=topoid)


        !----------------
        ! Define some fields. Vector and scalar fields, with different
        ! dimensions because the interpolation routines are hard-coded for some
        ! and we want to test them all!
        !----------------
        CALL VField1%create(2,info,name='vecField1')
        CALL VField1%discretize_on(Part1,info)
        CALL VField2%create(3,info,name='vecField2')
        CALL VField2%discretize_on(Part1,info)
        CALL VField3%create(4,info,name='vecField3')
        CALL VField3%discretize_on(Part1,info)
        CALL VField4%create(5,info,name='vecField4')
        CALL VField4%discretize_on(Part1,info)

        CALL SField1%create(1,info,name='scaField1')
        CALL SField1%discretize_on(Part1,info)
        CALL SField2%create(1,info,name='scaField2')
        CALL SField2%discretize_on(Part1,info)
        CALL SField3%create(1,info,name='scaField3')
        CALL SField3%discretize_on(Part1,info)
        CALL Vol%create(1,info,name='Part_Volume')
        CALL Vol%discretize_on(Part1,info)

        !----------------
        ! Initialize the fields with test functions (polynomials of orders 0,1
        ! and 2), to test interpolants of orders up to 3.
        !----------------
        foreach p in particles(Part1) with positions(x) vec_fields(V1=VField1,V2=VField2,V3=VField3,V4=Vfield4) sca_fields(S1=SField1,S2=SField2,S3=SField3,Vol=Vol)
                Vol_p   = 1.0_MK / np_global
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
        CALL Part1%map_ghosts(info)

        end_subroutine()
    end setup
!----------------------------------------------


!--------------- teardown ---------------------
    teardown
        CALL Mesh1%destroy(info)
        CALL SField1%destroy(info)
        CALL SField2%destroy(info)
        CALL SField3%destroy(info)
        CALL VField1%destroy(info)
        CALL VField2%destroy(info)
        CALL VField3%destroy(info)
        CALL VField4%destroy(info)
        CALL Vol%destroy(info)
        CALL Part1%destroy(info)

        NULLIFY(topo)
    end teardown
!----------------------------------------------

    test part_remesh

        start_subroutine("remesh")

        !CALL ppm_vtk_particles("output_before",Part1,info)
        !Assert_Equal(info,0)

        !----------------
        ! Remesh the particles
        ! (this performs the p2m interpolation as well)
        !----------------
        CALL Part1%interp_to_mesh(Mesh1,VField1,ppm_param_rmsh_kernel_mp4,info)
            Assert_Equal(info,0)

        CALL Part1%remesh(Mesh1,ppm_param_rmsh_kernel_mp4,info)
            Assert_Equal(info,0)

        !CALL ppm_vtk_particles("output_after",Part1,info)
        !Assert_Equal(info,0)

        tol = 1e-2

        !----------------
        ! Define a cutoff distance from the sides of the patch. Particles that
        ! are in the patch but too close to the sides will not receive anything
        ! from the patch during interpolation (except if there are some periodic
        ! boundaries, of course...)
        !----------------
        cutoff = REAL(Mesh1%ghostsize(1:ndim),ppm_kind_double)* Mesh1%h(1:ndim)

        voln = 1.0_MK / PRODUCT(Mesh1%Nm(1:ndim)-1)

        !----------------
        !Loop through all the mesh nodes and check that the values of the field
        !have been interpolated EXACTLY,
        !as should be the case for constant, linear and quadratic functions
        !with a 3rd order interpolating kernel (like Mp4)
        !----------------
        foreach p in particles(Part1) with positions(x) vec_fields(V1=VField1,V2=VField2,V3=VField3,V4=Vfield4) sca_fields(S1=SField1,S2=SField2,S3=SField3,Vol=Vol)
                Vol_p   = 1.0_MK / np_global
                Assert_Equal_Within(S1_p, f_cst(x_p(1:ndim),ndim)* Vol_p,tol)
                Assert_Equal_Within(S2_p, f_lin(x_p(1:ndim),ndim)* Vol_p,tol)
                Assert_Equal_Within(S3_p, f_sq(x_p(1:ndim),ndim) * Vol_p,tol)

                Assert_Equal_Within(V1_p(1), f_cst(x_p(1:ndim),ndim)* Vol_p,tol)
                Assert_Equal_Within(V1_p(2), f_sq(x_p(1:ndim),ndim) * Vol_p,tol)

                Assert_Equal_Within(V2_p(1), f_cst(x_p(1:ndim),ndim)* Vol_p,tol)
                Assert_Equal_Within(V2_p(2), f_lin(x_p(1:ndim),ndim)* Vol_p,tol)
                Assert_Equal_Within(V2_p(3), f_sq(x_p(1:ndim),ndim) * Vol_p,tol)

                Assert_Equal_Within(V3_p(1), f_cst(x_p(1:ndim),ndim)* Vol_p,tol)
                Assert_Equal_Within(V3_p(2), f_lin(x_p(1:ndim),ndim)* Vol_p,tol)
                Assert_Equal_Within(V3_p(3), f_sq(x_p(1:ndim),ndim) * Vol_p,tol)
                Assert_Equal_Within(V3_p(4), f_sq(x_p(1:ndim),ndim) * Vol_p,tol)

                Assert_Equal_Within(V4_p(1), f_cst(x_p(1:ndim),ndim)* Vol_p,tol)
                Assert_Equal_Within(V4_p(2), f_lin(x_p(1:ndim),ndim)* Vol_p,tol)
                Assert_Equal_Within(V4_p(3), f_sq(x_p(1:ndim),ndim) * Vol_p,tol)
                Assert_Equal_Within(V4_p(4), f_sq(x_p(1:ndim),ndim) * Vol_p,tol)
                Assert_Equal_Within(V4_p(5), f_sq(x_p(1:ndim),ndim) * Vol_p,tol)
        end foreach
#ifdef __MPI
        CALL MPI_BARRIER(comm,info)
#endif
        end_subroutine()
        !check that we are leaving the test without error
        Assert_Equal(info,0)
    end test

    test part_remesh_split

        start_subroutine("remesh_split")

        !CALL ppm_vtk_particles("output_before",Part1,info)
        !Assert_Equal(info,0)

        !----------------
        ! Remesh the particles
        ! (this performs the p2m interpolation as well)
        !----------------
        CALL Part1%interp_to_mesh_all(Mesh1,ppm_param_rmsh_kernel_mp4,info)
            Assert_Equal(info,0)

        !here, we could do something, like finite differences on the mesh
        ! (TODO: add a test for that)

        CALL Part1%recreate_from_mesh(Mesh1,info)
            Assert_Equal(info,0)

        !CALL ppm_vtk_particles("output_after",Part1,info)
        !Assert_Equal(info,0)

        tol = 1e-2

        !----------------
        ! Define a cutoff distance from the sides of the patch. Particles that
        ! are in the patch but too close to the sides will not receive anything
        ! from the patch during interpolation (except if there are some periodic
        ! boundaries, of course...)
        !----------------
        cutoff = REAL(Mesh1%ghostsize(1:ndim),ppm_kind_double)* Mesh1%h(1:ndim)

        voln = 1.0_MK / PRODUCT(Mesh1%Nm(1:ndim)-1)

        !----------------
        !Loop through all the mesh nodes and check that the values of the field
        !have been interpolated EXACTLY,
        !as should be the case for constant, linear and quadratic functions
        !with a 3rd order interpolating kernel (like Mp4)
        !----------------
        foreach p in particles(Part1) with positions(x) vec_fields(V1=VField1,V2=VField2,V3=VField3,V4=Vfield4) sca_fields(S1=SField1,S2=SField2,S3=SField3,Vol=Vol)
                Vol_p   = 1.0_MK / np_global
                Assert_Equal_Within(S1_p, f_cst(x_p(1:ndim),ndim)* Vol_p,tol)
                Assert_Equal_Within(S2_p, f_lin(x_p(1:ndim),ndim)* Vol_p,tol)
                Assert_Equal_Within(S3_p, f_sq(x_p(1:ndim),ndim) * Vol_p,tol)

                Assert_Equal_Within(V1_p(1), f_cst(x_p(1:ndim),ndim)* Vol_p,tol)
                Assert_Equal_Within(V1_p(2), f_sq(x_p(1:ndim),ndim) * Vol_p,tol)

                Assert_Equal_Within(V2_p(1), f_cst(x_p(1:ndim),ndim)* Vol_p,tol)
                Assert_Equal_Within(V2_p(2), f_lin(x_p(1:ndim),ndim)* Vol_p,tol)
                Assert_Equal_Within(V2_p(3), f_sq(x_p(1:ndim),ndim) * Vol_p,tol)

                Assert_Equal_Within(V3_p(1), f_cst(x_p(1:ndim),ndim)* Vol_p,tol)
                Assert_Equal_Within(V3_p(2), f_lin(x_p(1:ndim),ndim)* Vol_p,tol)
                Assert_Equal_Within(V3_p(3), f_sq(x_p(1:ndim),ndim) * Vol_p,tol)
                Assert_Equal_Within(V3_p(4), f_sq(x_p(1:ndim),ndim) * Vol_p,tol)

                Assert_Equal_Within(V4_p(1), f_cst(x_p(1:ndim),ndim)* Vol_p,tol)
                Assert_Equal_Within(V4_p(2), f_lin(x_p(1:ndim),ndim)* Vol_p,tol)
                Assert_Equal_Within(V4_p(3), f_sq(x_p(1:ndim),ndim) * Vol_p,tol)
                Assert_Equal_Within(V4_p(4), f_sq(x_p(1:ndim),ndim) * Vol_p,tol)
                Assert_Equal_Within(V4_p(5), f_sq(x_p(1:ndim),ndim) * Vol_p,tol)
        end foreach

#ifdef __MPI
        CALL MPI_BARRIER(comm,info)
#endif

        end_subroutine()
        !check that we are leaving the test without error
        Assert_Equal(info,0)
    end test

!-------------------------------------------------------------
! test function
!-------------------------------------------------------------
PURE FUNCTION f_cst(pos,ndim) RESULT(res)
    REAL(MK)                              :: res
    INTEGER                 ,  INTENT(IN) :: ndim
    REAL(MK), DIMENSION(ndim), INTENT(IN) :: pos

    res =  42.0_MK
END FUNCTION

PURE FUNCTION f_lin(pos,ndim) RESULT(res)
    REAL(MK)                              :: res
    INTEGER                 ,  INTENT(IN) :: ndim
    REAL(MK), DIMENSION(ndim), INTENT(IN) :: pos

    res =  1.337_MK + pos(1) + 10._MK*pos(2) + 100._MK*pos(ndim)
END FUNCTION

PURE FUNCTION f_sq(pos,ndim) RESULT(res)
    REAL(MK)                              :: res
    INTEGER                 ,  INTENT(IN) :: ndim
    REAL(MK), DIMENSION(ndim), INTENT(IN) :: pos

    res =  pos(1)**2 + 10._MK*pos(2)**2 + 100._MK*pos(ndim)**2
END FUNCTION

!!! check whether a particle is within a patch and more than a cutoff
!!! distance away from its boundaries.
PURE FUNCTION is_well_within(pos,patch,cutoff,ndim) RESULT(res)
    LOGICAL                               :: res
    REAL(MK), DIMENSION(ndim), INTENT(IN) :: pos
    REAL(MK), DIMENSION(2*ndim),INTENT(IN):: patch
    REAL(MK), DIMENSION(ndim), INTENT(IN) :: cutoff
    INTEGER                 ,  INTENT(IN) :: ndim

    res = ALL(pos(1:ndim).GE.(patch(1:ndim)+cutoff(1:ndim)))
    res = res .AND. ALL(pos(1:ndim).LE.(patch(ndim+1:2*ndim)-cutoff(1:ndim)))

END FUNCTION

end test_suite
