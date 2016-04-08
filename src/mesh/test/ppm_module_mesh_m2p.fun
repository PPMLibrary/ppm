test_suite ppm_module_mesh_m2p

USE ppm_module_mesh_typedef
USE ppm_module_topo_typedef
USE ppm_module_field_typedef
USE ppm_module_particles_typedef
USE ppm_module_mktopo

INTEGER, PARAMETER              :: debug = 0
INTEGER, PARAMETER              :: mk = KIND(1.0d0) !KIND(1.0e0)
REAL(MK),PARAMETER              :: pi = ACOS(-1._MK)
INTEGER,PARAMETER               :: ndim=3
INTEGER                         :: decomp,assig,tolexp
INTEGER                         :: info,comm,rank,nproc
REAL(MK)                        :: tol
INTEGER                         :: topoid=-1
REAL(MK),DIMENSION(:  ),POINTER :: min_phys => NULL()
REAL(MK),DIMENSION(:  ),POINTER :: max_phys => NULL()

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

!---------------- init -----------------------

    init

        USE ppm_module_topo_typedef
        USE ppm_module_init

        ALLOCATE(min_phys(ndim),max_phys(ndim),&
        & ighostsize(ndim),nm(ndim),h(ndim))

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
        CALL ppm_init(ndim,MK,tolexp,0,debug,info,99)

        CALL RANDOM_SEED(size=seedsize)
        ALLOCATE(seed(seedsize))
        DO i=1,seedsize
           seed(i)=9+i*(rank+1)
        ENDDO
        CALL random_seed(put=seed)

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
    end setup
!----------------------------------------------

!--------------- teardown ---------------------
    teardown
        NULLIFY(topo)
    end teardown
!----------------------------------------------

    test mesh_to_part_interp
        TYPE(ppm_t_field) :: VField1,VField2,VField3,VField4
        TYPE(ppm_t_field) :: SField1,SField2,SField3
        TYPE(ppm_t_particles_d),TARGET :: Part1

        REAL(ppm_kind_double),DIMENSION(ndim) :: x0

        REAL(ppm_kind_double),DIMENSION(ndim) :: pos
        REAL(ppm_kind_double),DIMENSION(ndim) :: cutoff
        INTEGER :: np_global = 500000

        start_subroutine("test_interp")

        !----------------
        ! make topology
        !----------------
        decomp = ppm_param_decomp_cuboid
        assig  = ppm_param_assign_internal
        topoid = 0
        sca_ghostsize = 0.07_MK
        CALL ppm_mktopo(topoid,decomp,assig,min_phys,max_phys,    &
        &    bcdef,sca_ghostsize,cost,info)
        Assert_Equal(info,0)


        !----------------
        ! Create Mesh
        !----------------
        Nm(1) = 35
        Nm(ndim) = 24
        Nm(2) = 65
        CALL Mesh1%create(topoid,offset,info,Nm=Nm, &
        &    ghostsize=ighostsize,name='Test_Mesh_1')
        Assert_Equal(info,0)

        !----------------
        ! Add a patch
        !----------------
        IF (ndim.EQ.2) THEN
           my_patch(1:4) = (/0.15_MK,0.10_MK,0.99_MK,0.7_MK/)
        ELSE
           my_patch(1:6) = (/0.15_MK,0.10_MK,0.5_MK,0.89_MK,0.7_MK,0.78_MK/)
        ENDIF
        CALL Mesh1%def_patch(my_patch,info)
        Assert_Equal(info,0)

        !----------------
        ! Create particles, from a grid + small random displacement
        !----------------
        CALL Part1%initialize(np_global,info,topoid=topoid,name="Part1")
        Assert_Equal(info,0)

        ALLOCATE(wp_2r(ndim,Part1%Npart))
        CALL RANDOM_NUMBER(wp_2r)
        wp_2r = (wp_2r - 0.5_MK) * Part1%h_avg * 1.15_MK

        CALL Part1%move(wp_2r,info)
        Assert_Equal(info,0)

        DEALLOCATE(wp_2r)

        !----------------
        ! Put particles back into the domain and global map them
        !----------------
        CALL Part1%apply_bc(info)
        Assert_Equal(info,0)

        CALL Part1%map(info,global=.TRUE.,topoid=topoid)
        Assert_Equal(info,0)


        !----------------
        ! Define some fields. Vector and scalar fields, with different
        ! DIMENSIONs because the interpolation routines are hard-coded for some
        ! and we want to test them all!
        !----------------
        CALL VField1%create(2,info,name='vecField1')
        CALL VField1%discretize_on(Mesh1,info)
        CALL VField2%create(3,info,name='vecField2')
        CALL VField2%discretize_on(Mesh1,info)
        CALL VField3%create(4,info,name='vecField3')
        CALL VField3%discretize_on(Mesh1,info)
        CALL VField4%create(5,info,name='vecField4')
        CALL VField4%discretize_on(Mesh1,info)

        CALL SField1%create(1,info,name='scaField1')
        CALL SField1%discretize_on(Mesh1,info)
        CALL SField2%create(1,info,name='scaField2')
        CALL SField2%discretize_on(Mesh1,info)
        CALL SField3%create(1,info,name='scaField3')
        CALL SField3%discretize_on(Mesh1,info)

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
        CALL Mesh1%interp_to_part(Part1,VField1,ppm_param_rmsh_kernel_mp4,info)
        Assert_Equal(info,0)
        CALL Mesh1%interp_to_part(Part1,VField2,ppm_param_rmsh_kernel_mp4,info)
        Assert_Equal(info,0)
        CALL Mesh1%interp_to_part(Part1,VField3,ppm_param_rmsh_kernel_mp4,info)
        Assert_Equal(info,0)
        CALL Mesh1%interp_to_part(Part1,VField4,ppm_param_rmsh_kernel_mp4,info)
        Assert_Equal(info,0)
        CALL Mesh1%interp_to_part(Part1,SField1,ppm_param_rmsh_kernel_mp4,info)
        Assert_Equal(info,0)
        CALL Mesh1%interp_to_part(Part1,SField2,ppm_param_rmsh_kernel_mp4,info)
        Assert_Equal(info,0)
        CALL Mesh1%interp_to_part(Part1,SField3,ppm_param_rmsh_kernel_mp4,info)
        Assert_Equal(info,0)

        tol = 1E-12

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
            IF (is_well_within(x_p(1:ndim),my_patch(1:2*ndim),cutoff,ndim)) THEN

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

            ENDIF
        end foreach
#ifdef __MPI
        CALL MPI_BARRIER(comm,info)
#endif

        CALL Mesh1%destroy(info)
        CALL SField1%destroy(info)
        CALL SField2%destroy(info)
        CALL SField3%destroy(info)
        CALL VField1%destroy(info)
        CALL VField2%destroy(info)
        CALL VField3%destroy(info)
        CALL VField4%destroy(info)
        CALL Part1%destroy(info)

        end_subroutine()
    end test

    !-------------------------------------------------------------
    ! test function
    !-------------------------------------------------------------
    PURE FUNCTION test_constant(pos,ndim) RESULT(res)
      IMPLICIT NONE
      INTEGER,                   INTENT(IN   ) :: ndim
      REAL(MK), DIMENSION(ndim), INTENT(IN   ) :: pos
      REAL(MK)                                 :: res

      res =  42.17_MK
    END FUNCTION

    PURE FUNCTION test_linear(pos,ndim) RESULT(res)
      IMPLICIT NONE
      INTEGER                 ,  INTENT(IN   ) :: ndim
      REAL(MK), DIMENSION(ndim), INTENT(IN   ) :: pos
      REAL(MK)                                 :: res
      res =  pos(1) + 10._MK*pos(2) + 100._MK*pos(ndim)
    END FUNCTION

    PURE FUNCTION test_quadratic(pos,ndim) RESULT(res)
      IMPLICIT NONE
      INTEGER                 ,  INTENT(IN   ) :: ndim
      REAL(MK), DIMENSION(ndim), INTENT(IN   ) :: pos
      REAL(MK)                                 :: res

      res =  pos(1)**2 + 10._MK*pos(2)**2 + 100._MK*pos(ndim)**2
    END FUNCTION

    !!! check whether a particle is within a patch and more than a cutoff
    !!! distance away from its boundaries.
    PURE FUNCTION is_well_within(pos,patch,cutoff,ndim) RESULT(res)
      IMPLICIT NONE
      REAL(MK), DIMENSION(ndim), INTENT(IN   ) :: pos
      REAL(MK), DIMENSION(2*ndim),INTENT(IN   ):: patch
      REAL(MK), DIMENSION(ndim), INTENT(IN   ) :: cutoff
      INTEGER                 ,  INTENT(IN   ) :: ndim
      LOGICAL                                  :: res

      res = ALL(pos(1:ndim).GE.(patch(1:ndim)+cutoff(1:ndim)))
      res = res .AND. ALL(pos(1:ndim).LE.(patch(ndim+1:2*ndim)-cutoff(1:ndim)))
    END FUNCTION

end test_suite
