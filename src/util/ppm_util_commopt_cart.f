      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_util_commopt_cart
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

      SUBROUTINE ppm_util_commopt_cart(topoid,info)
      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_topo_typedef
      USE ppm_module_substart
      USE ppm_module_write
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_mpi
      IMPLICIT NONE

      INTEGER, PARAMETER :: MK = ppm_kind_double
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER, INTENT(IN   ) :: topoid
      !!! The topology to be optimized
      INTEGER, INTENT(INOUT) :: info
      !!! return status, 0 on success

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      TYPE(ppm_t_topo), POINTER :: topo

      REAL(MK)                  :: t0
      REAL(MK), DIMENSION(3,26) :: displ
      REAL(MK), DIMENSION(3,3)  :: rotxy,rotyz
      REAL(MK), DIMENSION(3)    :: rdispl,dsplc
      REAL(MK), DIMENSION(4)    :: angles

      INTEGER, DIMENSION(3,8)    :: lcoords
      INTEGER, DIMENSION(3)      :: coords, lcoords2
      INTEGER, DIMENSION(3,26,8) :: commsequence
      INTEGER, DIMENSION(2)      :: ldu
      INTEGER                    :: idir,icpu,iangle,jangle,lrank
      ! processor ranks
      INTEGER                    :: idirect,lidir
      INTEGER                    :: iopt
      INTEGER                    :: i,j,k
#ifdef __MPI
      INTEGER, DIMENSION(3)      :: ndims
#endif

      CHARACTER(LEN=ppm_char) :: caller="ppm_util_commopt_cart"

      LOGICAL, DIMENSION(:), ALLOCATABLE :: ichecked
#ifdef __MPI
      LOGICAL, DIMENSION(3)              :: periods
#endif
      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
         CALL check
         IF (info .NE. 0) GOTO 9999
      ENDIF

      topo => ppm_topo(topoid)%t
      !-------------------------------------------------------------------------
      !  Check if there are more than 1 processor. If not, we are done
      !-------------------------------------------------------------------------
      IF (ppm_nproc.LT.2) THEN
         !---------------------------------------------------------------------
         !  Allocate memory for communication protocols
         !---------------------------------------------------------------------
         iopt  =ppm_param_alloc_grow_preserve
         ldu(1)=1
         CALL ppm_alloc(topo%icommseq,ldu,iopt,info)
         or_fail_alloc("communication sequence PPM_ICOMMSEQ",ppm_error=ppm_error_fatal)

         !---------------------------------------------------------------------
         !  Set the trivial protocol: only myself
         !---------------------------------------------------------------------
         topo%ncommseq = 1
         topo%icommseq(1) = ppm_rank
         GOTO 9999
      ENDIF

#ifdef __MPI

      !-----------------------------------------------------
      ! Retrieves Cartesian topology information associated with a communicator
      !
      ! Get our coordinates and the number of cpus per dimension
      ! (need special care if one of the dim. is 1)
      !-----------------------------------------------------
      CALL MPI_Cart_get(ppm_comm,ppm_dim,ndims,periods,coords,info)
      or_fail_MPI("MPI_CART_GET")

      IF (ppm_debug .GT. 2) THEN
         IF (ppm_dim.EQ.2) THEN
            stdout_f('(A,2I4,A,2I4,A)',"get neighbors for MPI-coordinates ",coords," (ndims = ",ndims,")")
         ELSE
            stdout_f('(A,3I4,A,3I4,A)',"get neighbors for MPI-coordinates ",coords," (ndims = ",ndims,")")
         ENDIF
      ENDIF

      !-----------------------------------------------------
      !  Allocate the checked array
      !-----------------------------------------------------
      ALLOCATE(ichecked(ppm_nproc),STAT=info)
      or_fail_alloc("ichecked array",ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      icpu = 1
      DO k=1,2
         DO j=1,2
            DO i=1,2
               lcoords(:,icpu) = (/i-1,j-1,k-1/)
               icpu = icpu + 1
            ENDDO
         ENDDO
      ENDDO
      rdispl    = (/1.0,0.0,0.0/)
      angles(1) = 0.0
      angles(2) = 0.25_MK * ppm_pi_d
      angles(3) = 0.5_MK  * ppm_pi_d
      angles(4) = 0.75_MK * ppm_pi_d
      idir = 0
      DO jangle=1,4
         DO iangle=1,4
            rotxy(:,1) = (/ COS(angles(jangle)), SIN(angles(jangle)),0.0_MK/)
            rotxy(:,2) = (/-SIN(angles(jangle)), COS(angles(jangle)),0.0_MK/)
            rotxy(:,3) = (/ 0.0_MK, 0.0_MK, 1.0_MK /)
            rotyz(:,1) = (/ 1.0_MK, 0.0_MK, 0.0_MK /)
            rotyz(:,2) = (/0.0_MK, COS(angles(iangle)),SIN(angles(iangle))/)
            rotyz(:,3) = (/0.0_MK,-SIN(angles(iangle)),COS(angles(iangle))/)
            dsplc = MATMUL(MATMUL(rotyz,rotxy),rdispl)
            idir = idir + 1
            IF (idir.LT.4) CYCLE
            displ(:,2*(idir-3)-1) = dsplc
            displ(:,2*(idir-3)  ) = -dsplc
         ENDDO
      ENDDO
      !-----------------------------------------------------
      !  truncate
      DO idir=1,26
         DO i=1,3
            IF (ABS(displ(i,idir)).LT.0.05_MK) THEN
               displ(i,idir) = 0.0_MK
            ENDIF
         ENDDO
      ENDDO
      !-----------------------------------------------------
      !  modulate
      DO icpu=1,8
         DO idir=1,26
            idirect = 1
            IF (ndims(1).GT.1 .AND. ABS(displ(1,idir)).GT.0.0) THEN
               idirect = 1
            ENDIF
            IF (ndims(2).GT.1 .AND. ABS(displ(2,idir)).GT.0.0) THEN
               idirect = 2
            ENDIF
            IF (ndims(3).GT.1 .AND. ABS(displ(3,idir)).GT.0.0) THEN
               idirect = 3
            ENDIF
            IF (MOD(lcoords(idirect,icpu),2).EQ.1) THEN
               IF (displ(1,idir).LT.0.0) THEN
                  commsequence(1,idir,icpu) = -1
               ELSE IF (displ(1,idir).GT.0.0) THEN
                  commsequence(1,idir,icpu) =  1
               ELSE
                  commsequence(1,idir,icpu) =  0
               ENDIF
               IF (displ(2,idir).LT.0.0) THEN
                  commsequence(2,idir,icpu) = -1
               ELSE IF (displ(2,idir).GT.0.0) THEN
                  commsequence(2,idir,icpu) =  1
               ELSE
                  commsequence(2,idir,icpu) =  0
               ENDIF
               IF (displ(3,idir).LT.0.0) THEN
                  commsequence(3,idir,icpu) = -1
               ELSE IF (displ(3,idir).GT.0.0) THEN
                  commsequence(3,idir,icpu) =  1
               ELSE
                  commsequence(3,idir,icpu) =  0
               ENDIF
            ELSE
               IF (displ(1,idir).LT.0.0) THEN
                  commsequence(1,idir,icpu) =  1
               ELSE IF (displ(1,idir).GT.0.0) THEN
                  commsequence(1,idir,icpu) = -1
               ELSE
                  commsequence(1,idir,icpu) =  0
               ENDIF
               IF (displ(2,idir).LT.0.0) THEN
                  commsequence(2,idir,icpu) =  1
               ELSE IF (displ(2,idir).GT.0.0) THEN
                  commsequence(2,idir,icpu) = -1
               ELSE
                  commsequence(2,idir,icpu) =  0
               ENDIF
               IF (displ(3,idir).LT.0.0) THEN
                  commsequence(3,idir,icpu) =  1
               ELSE IF (displ(3,idir).GT.0.0) THEN
                  commsequence(3,idir,icpu) = -1
               ELSE
                  commsequence(3,idir,icpu) =  0
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      !-----------------------------------------------------
      ! Determines process coords in Cartesian topology given rank in group.
      !
      ! we have got the basic pattern now, now evaluate it for this rank
      !-----------------------------------------------------
      CALL MPI_Cart_coords(ppm_comm,ppm_rank,ppm_dim,coords,info)
      or_fail_MPI("MPI_Cart_coords")

      ichecked = .FALSE.

      topo%ncommseq = MIN(27,ppm_nproc)
      iopt   = ppm_param_alloc_grow_preserve
      ldu(1) = topo%ncommseq
      CALL ppm_alloc(topo%icommseq,ldu,iopt,info)
      or_fail_alloc("final communication sequence PPM_ICOMMSEQ",ppm_error=ppm_error_fatal)

      !-----------------------------------------------------
      !  what is my icpu
      icpu = 1*MOD(coords(1),2)+2*MOD(coords(2),2)+4*MOD(coords(3),2)+1
      topo%icommseq(1) = ppm_rank
      ichecked(ppm_rank+1) = .TRUE.
      lidir = 1

      DO idir=1,26
         lcoords2 = coords + commsequence(:,idir,icpu)
         ! stdout(coords,"/",lcoords2)
         lrank = -1

         !  Determines process rank in communicator given Cartesian location.
         CALL MPI_Cart_rank(ppm_comm,lcoords2,lrank,info)
         ! stdout("rank =>",lrank,"(",info,")")

         IF (ichecked(lrank+1)) CYCLE

         topo%icommseq(lidir+1) = lrank
         ichecked(lrank+1) = .TRUE.
         lidir = lidir + 1
      ENDDO

      topo%ncommseq = lidir

      !IF (ppm_rank.EQ.0) THEN
      !   stdout("ppm_ncommseq is ",'ppm_ncommseq(topoid)',lidir)
      !   DO idir=1,ppm_ncommseq(topoid)
      !      stdout(" entry ",'ppm_icommseq(idir,topoid)')
      !   ENDDO
      !ENDIF

      !  done
      !-----------------------------------------------------
      !-------------------------------------------------------------------------
      !  Everybody gets the memory needed
      !-------------------------------------------------------------------------
      !stdout("ppm_icommseq ",'ppm_icommseq(1:9,topoid)')
      !stdout("ppm_icommseq ",'ppm_icommseq(10:18,topoid)')
      !stdout("ppm_icommseq ",'ppm_icommseq(19:27,topoid)')

      !-------------------------------------------------------------------------
      !  Mark this topology as done
      !-------------------------------------------------------------------------
      topo%isoptimized = .TRUE.

      !-------------------------------------------------------------------------
      !  Deallocate temporary storage
      !-------------------------------------------------------------------------
      DEALLOCATE(ichecked,STAT=info)
      or_fail_dealloc("number of comm.rounds PPM_NCOMMSEQ",ppm_error=ppm_error_fatal)
#endif

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
    9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
    CONTAINS
      SUBROUTINE check
        IMPLICIT NONE
        LOGICAL :: valid
        CALL ppm_check_topoid(topoid,valid,info)
        IF (.NOT. valid) THEN
           fail("topoid out of range",exit_point=8888)
        ENDIF
    8888 CONTINUE
      END SUBROUTINE check
    END SUBROUTINE ppm_util_commopt_cart
