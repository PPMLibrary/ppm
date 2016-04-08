      !-------------------------------------------------------------------------
      !     Test Case   :                   ppm_test_interp_p2m
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
program ppm_test_interp_p2m

    USE ppm_module_typedef
    USE ppm_module_mktopo
    USE ppm_module_topo_get
    USE ppm_module_interp_m2p
    USE ppm_module_interp_p2m
    USE ppm_module_init
    USE ppm_module_finalize
    USE ppm_module_map
    USE ppm_module_mpi
    IMPLICIT NONE

#include "../../ppm_define.h"

    INTEGER, PARAMETER              :: debug = 0
    INTEGER, PARAMETER              :: MK = ppm_kind_double
    INTEGER                         :: ndim,nspec
    INTEGER                         :: decomp
    INTEGER                         :: assig
    INTEGER                         :: tolexp
    REAL(MK)                        :: tol
    INTEGER                         :: info
    INTEGER                         :: topoid,meshid
    REAL(MK),DIMENSION(:,:),POINTER :: xp,wp
    REAL(MK),DIMENSION(:  ),POINTER :: min_phys,max_phys,h
    INTEGER, DIMENSION(:  ),POINTER :: ghostsize
    INTEGER                         :: i,j,p_i,ai,aj
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
    !----------------
    ! setup
    !----------------
    tol = 10.0_MK*EPSILON(1.0_MK)
    tolexp = INT(LOG10(EPSILON(1.0_MK)))
    ndim = 2
    nspec = 1
    data ((alpha(ai,aj), ai=1,2), aj=1,nmom) /0,0, 1,0, 0,1, 2,0, 0,2, &
   &                                          1,1, 3,0, 0,3, 2,1, 1,2/

    ALLOCATE(min_phys(ndim),max_phys(ndim),ghostsize(ndim),&
   &         nm(ndim),h(ndim),field_x(ndim),STAT=info)

    DO i=1,ndim
        min_phys(i) = 0.0_MK
        max_phys(i) = 1.0_MK
        ghostsize(i) = 2
    ENDDO
    bcdef(1:6) = ppm_param_bcdef_freespace

    NULLIFY(xp)
    NULLIFY(wp)

#ifdef __MPI
    CALL MPI_Init(info)
#endif

    CALL ppm_init(ndim,MK,tolexp,0,debug,info,99)

    !----------------
    ! create particles
    !----------------
    np = 400
    mp = 0

    ALLOCATE(xp(ndim,np),wp(nspec,np),STAT=info)
    !CALL RANDOM_SEED(put=(/17,42/))
    CALL RANDOM_NUMBER(xp)
    wp = 0.0_MK


    !----------------
    ! make topology
    !----------------
    decomp = ppm_param_decomp_user_defined
    assig = ppm_param_assign_user_defined

    topoid = 0
    meshid = -1

    ALLOCATE(minsub(ndim,1),maxsub(ndim,1),sub2proc(1),cost(1),STAT=info)

    cost(1) = 1.0_MK

    DO i=1,ndim
        minsub(i,1)=min_phys(i)
        maxsub(i,1)=max_phys(i)
    ENDDO

    sub2proc(1) = 0

    ALLOCATE(nm(ndim),STAT=info)
    DO i=1,ndim
        nm(i) = 32
    ENDDO

    CALL ppm_mktopo(topoid,meshid,xp,np,decomp,assig,min_phys,max_phys,bcdef, &
    &               ghostsize,cost,istart,ndata,nm,info,user_minsub=minsub,   &
    &               user_maxsub=maxsub,user_nsubs=1,user_sub2proc=sub2proc)

    ALLOCATE(field_wp(nspec,(1-ghostsize(1)):(ndata(1,1)+ghostsize(1)),         &
    &        (1-ghostsize(2)):(ndata(2,1)+ghostsize(2)),1),STAT=info) ! 2d
    !ALLOCATE(field_wp(nspec,ndata(1,1),ndata(2,1),ndata(3,1),1),STAT=info) ! 3d

    DO i=1,ndim
        h(i) = (max_phys(i) - min_phys(i)) / REAL(ndata(i,1)-1,mk)
    ENDDO

    NULLIFY(topo)
    CALL ppm_topo_get(topoid,topo,info)
    print *,topoid,meshid

    print *, 'finished initialization'

    CALL ppm_map_part_global(topoid,xp,np,info) ! positions
    CALL ppm_map_part_push(wp,nspec,np,info)    ! strengths
    CALL ppm_map_part_send(np,mp,info)          ! send
    CALL ppm_map_part_pop(wp,nspec,np,mp,info)  ! strengths
    CALL ppm_map_part_pop(xp,ndim,np,mp,info)   ! positions
    IF (info .NE. 0) then
       print *, 'Failed to do global mapping'
       goto 8000
    ENDIF
    np = mp

    maxm3 = 0.0_MK
    do p_i=1,np
        !----------------
        ! p --> m
        !----------------
        wp(1,p_i) = 1.0_MK
        CALL ppm_interp_p2m(topoid,meshid,xp,np,wp,1,kernel,ghostsize,field_wp,info)

        !----------------
        ! test p --> m
        !----------------
        f_moments = 0.0_MK
        p_moments = 0.0_MK
        DO j = 1-ghostsize(2), ndata(2,1)+ghostsize(2)
            DO i = 1-ghostsize(1),ndata(1,1)+ghostsize(1)
                field_x(1) = min_phys(1) + h(1)*REAL(i-1,mk)
                field_x(2) = min_phys(2) + h(2)*REAL(j-1,mk)
                do aj = 1,nmom
                   f_moments(aj) = f_moments(aj) + field_wp(1,i,j,1)* &
    &                            field_x(1)**alpha(1,aj)*field_x(2)**alpha(2,aj)
                ENDDO
            ENDDO
        ENDDO
        do aj = 1,nmom
           p_moments(aj) = xp(1,p_i)**alpha(1,aj)*xp(2,p_i)**alpha(2,aj)
        ENDDO
        do aj = 1,6
            IF (ABS(f_moments(aj) - p_moments(aj)) .GT. tol) then
                print *, 'particle pos:',     xp(:,p_i)
                print *, 'failed at moment: ', aj
                print *, 'field moments: ',   f_moments
                print *, 'particle moments: ',p_moments
                stop 'ERROR: p2m interpolation: moments not conserved.'
            ENDIF
        ENDDO
        do aj = 7,10
            IF (ABS(f_moments(aj) - p_moments(aj)) .GT. maxm3) then
                maxm3 = ABS(f_moments(aj) - p_moments(aj))
            ENDIF
        ENDDO

        wp(1,p_i) = 0.0_MK
    ENDDO

    print *, 'Maximum 3rd moment difference / h^3', maxm3/h**3
    !----------------
    ! m --> p
    !----------------


    !----------------
    ! teardown...
    !----------------
    print *, 'cleaning up'
8000 CONTINUE

    CALL ppm_finalize(info)

#ifdef __MPI
    CALL MPI_Finalize(info)
#endif

    DEALLOCATE(xp,wp,field_wp,min_phys,max_phys,ghostsize,nm)

    print *, 'done.'

END PROGRAM ppm_test_interp_p2m
