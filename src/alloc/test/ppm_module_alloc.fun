test_suite ppm_module_alloc



#ifdef __MPI
    INCLUDE "mpif.h"
#endif

integer, parameter              :: debug = 0
integer, parameter              :: mk = kind(1.0d0) !kind(1.0e0)
integer,parameter               :: ndim=2
integer                         :: tolexp
integer                         :: nproc,rank,comm
real(mk),dimension(:,:),pointer :: xp => NULL()
real(mk),dimension(:),pointer   :: wp => NULL()
real(mk)                        :: t0
real(mk)                        :: tol
integer                         :: info

    init

        use ppm_module_init

        tol = epsilon(1.0_MK)
        tolexp = int(log10(epsilon(1.0_MK)))

        nullify(xp)

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


    finalize
        use ppm_module_finalize

        call ppm_finalize(info)

    end finalize



    test alloc
    ! test simple allocation
    use ppm_module_data
    integer, dimension(1) :: lda
    integer               :: iopt
    integer               :: info

    iopt = ppm_param_alloc_fit
    lda(1) = 1024
    call ppm_alloc(wp,lda,iopt,info)
    assert_equal(info,0)

    end test

    test dealloc
    ! test simple allocation
    use ppm_module_data
    integer, dimension(1) :: lda
    integer               :: iopt
    integer               :: info

    iopt = ppm_param_dealloc
    lda(1) = 1024
    call ppm_alloc(wp,lda,iopt,info)
    assert_equal(info,0)

    end test



end test_suite
