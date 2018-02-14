      !-------------------------------------------------------------------------
      !     Test Case   :                   ppm_test_CASE
      !-------------------------------------------------------------------------
      ! Copyright (c) 2012 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich),
      !                    Center for Fluid Dynamics (DTU)
      !
      !
      ! This file is part of the Parallel Particle Mesh Library (PPM).
      !
      ! PPM is free software: you can redistribute it and/or modify
      ! it under the terms of the GNU Lesser General Public License
      ! as published by the Free Software Foundation, either
      ! version 3 of the License, or (at your option) any later
      ! version.
      !
      ! PPM is distributed in the hope that it will be useful,
      ! but WITHOUT ANY WARRANTY; without even the implied warranty of
      ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
      ! GNU General Public License for more details.
      !
      ! You should have received a copy of the GNU General Public License
      ! and the GNU Lesser General Public License along with PPM. If not,
      ! see <http://www.gnu.org/licenses/>.
      !
      ! Parallel Particle Mesh Library (PPM)
      ! ETH Zurich
      ! CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
program ppm_test_CASE

    USE ppm_module_core
    USE ppm_module_typedef
    USE ppm_module_mktopo
    USE ppm_module_init
    USE ppm_module_finalize
    USE ppm_module_mpi
    IMPLICIT NONE

#include "../../ppm_define.h"

    INTEGER, PARAMETER              :: debug = 0
    INTEGER, PARAMETER              :: MK = ppm_kind_double
    REAL(MK),parameter              :: pi = 3.1415926535897931_MK
    INTEGER                         :: ndim,nspec
    INTEGER                         :: decomp
    INTEGER                         :: assig
    INTEGER                         :: tolexp
    REAL(MK)                        :: tol
    INTEGER                         :: info, comm
    INTEGER                         :: topoid,meshid
    INTEGER,  parameter             :: ngrid = 513
    INTEGER,  parameter             :: npgrid = 1025
    REAL(MK),DIMENSION(:,:),POINTER :: xp,wp
    REAL(MK),DIMENSION(:  ),POINTER :: min_phys,max_phys,h,p_h
    INTEGER, DIMENSION(:  ),POINTER :: ghostsize
    INTEGER                         :: i,j,ai,aj,p_i
    INTEGER, DIMENSION(6)           :: bcdef
    REAL(MK),DIMENSION(:  ),POINTER :: cost
    INTEGER, DIMENSION(:,:),POINTER :: istart,ndata
    INTEGER, DIMENSION(:  ),POINTER :: nm
    REAL(MK),DIMENSION(:,:),POINTER :: minsub,maxsub
    INTEGER, DIMENSION(:  ),POINTER :: sub2proc
    INTEGER                         :: np,mp
    INTEGER, PARAMETER              :: kernel = ppm_param_rmsh_kernel_mp4
    REAL(MK),DIMENSION(:,:,:,:  ), POINTER :: field_wp ! 2d  field_up(ldn,i,j,isub)
    !REAL(MK),DIMENSION(:,:,:,:,:), POINTER :: field_up ! 3d  field_up(ldn,i,j,k,isub)
    REAL(MK),DIMENSION(:  ),POINTER :: field_x
    TYPE(ppm_t_topo), POINTER       :: topo
    REAL(MK)                        :: maxm3
    INTEGER, PARAMETER              :: nmom = 10
    INTEGER, DIMENSION(2,nmom)      :: alpha
    REAL(MK),DIMENSION(nmom)        :: f_moments, p_moments
    REAL(MK),DIMENSION(2)           :: sigma,mu
    !================
    ! setup
    !================
    tol = 10.0_MK*EPSILON(1.0_MK)
    tolexp = INT(LOG10(EPSILON(1.0_MK)))
    ndim = 2
    nspec = 1

    ALLOCATE(min_phys(ndim),max_phys(ndim),ghostsize(ndim),&
   &         nm(ndim),h(ndim),p_h(ndim),field_x(ndim),STAT=info)

    DO i=1,ndim
        min_phys(i) = 0.0_MK
        max_phys(i) = 1.0_MK
        ghostsize(i) = 2
    ENDDO
    bcdef(1:6) = -1 ! BOUNDARY CONDITION

    NULLIFY(xp)
    NULLIFY(wp)

#ifdef __MPI
    CALL MPI_Init(info)
    comm = MPI_COMM_WORLD
#endif

    CALL ppm_init(ndim,MK,tolexp,comm,debug,info,99)

    !----------------
    ! create particles
    !----------------
    np = npgrid**2
    mp = 0

    ALLOCATE(xp(ndim,np),wp(nspec,np),STAT=info)
    xp = 0.0_MK
    wp = 0.0_MK


    ALLOCATE(nm(ndim),STAT=info)
    DO i=1,ndim
        nm(i) = ngrid
    ENDDO

    !----------------
    ! make topology
    !----------------

    topoid = 0
    meshid = -1

    CALL ppm_mktopo(topoid,meshid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
    &               ghostsize,cost,istart,ndata,nm,info)

    !----------------
    ! setup mesh data
    !----------------


    print *, 'finished initialization'
    !================
    ! body
    !================

    !---------------------------
    ! CALL routines to be tested
    !---------------------------

    !----------------
    ! test output
    !----------------




    !================
    ! teardown...
    !================
    print *, 'cleaning up'
8000 CONTINUE

    CALL ppm_finalize(info)

#ifdef __MPI
    CALL MPI_Finalize(info)
#endif

    DEALLOCATE(xp,wp,field_wp,min_phys,max_phys,ghostsize,nm)

    print *, 'done.'

END PROGRAM ppm_test_CASE
