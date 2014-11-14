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
      !!! Determine an approximately optimal communication
      !!! sequence for each processor to SendRecv data with
      !!! its neighbors only. Such that: no conflicts occur
      !!! (i.e. A wants to send to B, but B is currently busy
      !!! receiving from C) and that a minimum number of
      !!! communication rounds are needed. This is done by
      !!! using the Vizing approximation to the minimal edge
      !!! coloring problem of the processor topology graph.
      !!!
      !!! [TIP]
      !!! This routine should be used for cartesian decompositions
      !!!
      !!! [NOTE]
      !!! The current implementation relies on external routine written in C++
      !!! and its Fortran wrapper.
      !!!
      !!! .References
      !!! *******************************************************
      !!! - V.G. Vizing, On an estimate of the chromatic class
      !!! of a p-graph. Discret. Analiz. 3, 1964, pp.25-30.
      !!! (In Russian).
      !!! *******************************************************
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
      USE ppm_module_check_id
      IMPLICIT NONE

      INTEGER, PARAMETER :: MK = ppm_kind_double
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: topoid
      !!! The topology to be optimized
      INTEGER                 , INTENT(INOUT) :: info
      !!! return status, 0 on success

      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3,8)               :: lcoords
      LOGICAL, DIMENSION(:), POINTER        :: ichecked => NULL()
      INTEGER, DIMENSION(3)                 :: coords, lcoords2
      REAL(MK),    DIMENSION(3,26)          :: displ
      INTEGER, DIMENSION(3,26,8)            :: commsequence
      REAL(MK),    DIMENSION(3,3)           :: rotxy,rotyz
      REAL(MK),    DIMENSION(3)             :: rdispl,dsplc
      REAL(MK),    DIMENSION(4)             :: angles
      INTEGER                               :: idir,icpu,iangle,jangle,lrank
      INTEGER                               :: idirect,lidir
      !-----------------------------------------------------
      REAL(MK)                              :: t0
      INTEGER, DIMENSION(2)                 :: ldu
      INTEGER                               :: iopt
      ! number of neighborhood relations in total
!       INTEGER                               :: nlinks
      ! DISTINCT neighbor pairs, i.e. 1<->2 and 2<->1 is the same and only
      ! listed once. even indices are first points, odd ones second ones of
      ! the same edges.
!       INTEGER, DIMENSION(:  ) , POINTER     :: ilinks      => NULL()
      ! optimal edge coloring determined. sequence of triples (p1,p2,c),...
      ! with p1 and p2 being the 2 vertices of each edge and c its color.
!       INTEGER, DIMENSION(:  ) , POINTER     :: optres      => NULL()
      ! number of neighbors of all every CPU. index: MPI rank
!       INTEGER, DIMENSION(:  ) , POINTER     :: nneighprocs => NULL()
      ! all neighbors of all processors. 1st index: neighbor nr., 2nd:
      ! processor rank
!       INTEGER, DIMENSION(:,:) , POINTER     :: ineighprocs => NULL()
      INTEGER                               :: i,j,maxneigh,isize,ii,isin,k
      ! processor ranks
      INTEGER                               :: p1,p2,isb,jsb,ksb
      ! min and max of assigned colors
      INTEGER                               :: mincolor,maxcolor
      LOGICAL                               :: valid
#ifdef __MPI
      ! MPI comm status
      INTEGER, DIMENSION(MPI_STATUS_SIZE)   :: status
      LOGICAL, DIMENSION(3)                 :: periods
      INTEGER, DIMENSION(3)                 :: ndims
#endif
      TYPE(ppm_t_topo),      POINTER        :: topo

      CHARACTER(len=256) :: mesg
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_util_commopt_cart',t0,info)

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
      IF (ppm_nproc .LT. 2) THEN
          !---------------------------------------------------------------------
          !  Allocate memory for communication protocols
          !---------------------------------------------------------------------
          iopt   = ppm_param_alloc_grow_preserve
          ldu(1) = 1
          CALL ppm_alloc(topo%icommseq,ldu,iopt,info)
          IF (info .NE. 0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_alloc,'ppm_util_commopt_cart',     &
     &            'communication sequence PPM_ICOMMSEQ',__LINE__,info)
              GOTO 9999
          ENDIF


          !---------------------------------------------------------------------
          !  Set the trivial protocol: only myself
          !---------------------------------------------------------------------
          topo%ncommseq = 1
          topo%icommseq(1) = ppm_rank
          GOTO 9999
      END IF

#ifdef __MPI

      !-----------------------------------------------------
      !  Get our coordinates and the number of cpus per dimension
      !  (need special care if one of the dim. is 1)
      !-----------------------------------------------------
      CALL MPI_CART_GET(ppm_comm, 3, ndims, periods, coords, info)
      IF (ppm_debug .GT. 2) THEN
         WRITE (mesg, '(A,3I3,A,3I3,A)') 'get neighbors for MPI-coordinates ',&
    &           coords, ' (ndims = ', ndims, ')'
         CALL ppm_write(ppm_rank,'ppm_util_commopt_cart',mesg,info)
      ENDIF

      !-----------------------------------------------------
      !  Allocate the checked array
      !-----------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldu(1) = ppm_nproc
      CALL ppm_alloc(ichecked,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_util_commopt_cart',     &
              &            'ichecked array',__LINE__,info)
         GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      icpu = 1
      DO k=1,2
         DO j=1,2
            DO i=1,2
               lcoords(:,icpu) = (/i-1,j-1,k-1/)
               icpu = icpu + 1
            END DO
         END DO
      END DO
      rdispl    = (/1.0,0.0,0.0/)
      angles(1) = 0.0
      angles(2) = 0.25 * ppm_pi_d
      angles(3) = 0.5  * ppm_pi_d
      angles(4) = 0.75 * ppm_pi_d
      idir = 0
      DO jangle=1,4
         DO iangle=1,4
            rotxy(:,1) = (/ DCOS(angles(jangle)), DSIN(angles(jangle)),0.0_MK/)
            rotxy(:,2) = (/-DSIN(angles(jangle)), DCOS(angles(jangle)),0.0_MK/)
            rotxy(:,3) = (/ 0.0_MK, 0.0_MK, 1.0_MK /)
            rotyz(:,1) = (/ 1.0_MK, 0.0_MK, 0.0_MK /)
            rotyz(:,2) = (/0.0_MK, DCOS(angles(iangle)),DSIN(angles(iangle))/)
            rotyz(:,3) = (/0.0_MK,-DSIN(angles(iangle)),DCOS(angles(iangle))/)
            dsplc = MATMUL(MATMUL(rotyz,rotxy),rdispl)
            idir = idir + 1
            IF(idir.LT.4) CYCLE
            displ(:,2*(idir-3)-1) = dsplc
            displ(:,2*(idir-3)  ) = -dsplc
         END DO
      END DO
      !-----------------------------------------------------
      !  truncate
      DO idir=1,26
         DO i=1,3
            IF(ABS(displ(i,idir)).LT.0.05) THEN
               displ(i,idir) = 0.0
            END IF
         END DO
      END DO
      !-----------------------------------------------------
      !  modulate
      DO icpu=1,8
         DO idir=1,26
            idirect = 1
            IF(ndims(1).GT.1 .AND. ABS(displ(1,idir)).GT.0.0) THEN
               idirect = 1
            END IF
            IF(ndims(2).GT.1 .AND. ABS(displ(2,idir)).GT.0.0) THEN
               idirect = 2
            END IF
            IF(ndims(3).GT.1 .AND. ABS(displ(3,idir)).GT.0.0) THEN
               idirect = 3
            END IF
            IF(MOD(lcoords(idirect,icpu),2).EQ.1) THEN
               IF(displ(1,idir).LT.0.0) THEN
                  commsequence(1,idir,icpu) = -1
               ELSEIF(displ(1,idir).GT.0.0) THEN
                  commsequence(1,idir,icpu) =  1
               ELSE
                  commsequence(1,idir,icpu) =  0
               END IF
               IF(displ(2,idir).LT.0.0) THEN
                  commsequence(2,idir,icpu) = -1
               ELSEIF(displ(2,idir).GT.0.0) THEN
                  commsequence(2,idir,icpu) =  1
               ELSE
                  commsequence(2,idir,icpu) =  0
               END IF
               IF(displ(3,idir).LT.0.0) THEN
                  commsequence(3,idir,icpu) = -1
               ELSEIF(displ(3,idir).GT.0.0) THEN
                  commsequence(3,idir,icpu) =  1
               ELSE
                  commsequence(3,idir,icpu) =  0
               END IF

            ELSE
               IF(displ(1,idir).LT.0.0) THEN
                  commsequence(1,idir,icpu) =  1
               ELSEIF(displ(1,idir).GT.0.0) THEN
                  commsequence(1,idir,icpu) = -1
               ELSE
                  commsequence(1,idir,icpu) =  0
               END IF
               IF(displ(2,idir).LT.0.0) THEN
                  commsequence(2,idir,icpu) =  1
               ELSEIF(displ(2,idir).GT.0.0) THEN
                  commsequence(2,idir,icpu) = -1
               ELSE
                  commsequence(2,idir,icpu) =  0
               END IF
               IF(displ(3,idir).LT.0.0) THEN
                  commsequence(3,idir,icpu) =  1
               ELSEIF(displ(3,idir).GT.0.0) THEN
                  commsequence(3,idir,icpu) = -1
               ELSE
                  commsequence(3,idir,icpu) =  0
               END IF
            END IF
         END DO
      END DO
      !-----------------------------------------------------
      !   we have got the basic pattern now
      !   now evaluate it for this rank
      CALL MPI_CART_COORDS(ppm_comm,ppm_rank,3,coords,info)
      !-----------------------------------------------------
      ichecked = .FALSE.

      topo%ncommseq = MIN(27,ppm_nproc)
      iopt   = ppm_param_alloc_grow_preserve
      ldu(1) = topo%ncommseq
      CALL ppm_alloc(topo%icommseq,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_util_commopt_cart',     &
              &        'final communication sequence PPM_ICOMMSEQ',__LINE__,info)
         GOTO 9999
      ENDIF
      !-----------------------------------------------------
      !  what is my icpu
      icpu = 1*MOD(coords(1),2)+2*MOD(coords(2),2)+4*MOD(coords(3),2)+1
      topo%icommseq(1) = ppm_rank
      ichecked(ppm_rank+1) = .TRUE.
      lidir = 1
      DO idir=1,26
         lcoords2 = coords + commsequence(:,idir,icpu)
         !WRITE(mesg,*) coords,'/',lcoords2
         !CALL ppm_write(ppm_rank,'commopt',mesg,info)
         lrank = -1
         CALL MPI_CART_RANK(ppm_comm,lcoords2,lrank,info)
         !WRITE(mesg,*) 'rank =>',lrank,'(',info,')'
         !CALL ppm_write(ppm_rank,'commopt',mesg,info)
         IF(ichecked(lrank+1)) CYCLE
         topo%icommseq(lidir+1) = lrank
         ichecked(lrank+1) = .TRUE.
         lidir = lidir + 1
      END DO
      topo%ncommseq = lidir
      ! IF(ppm_rank.EQ.0) THEN
!          WRITE(mesg,*) 'ppm_ncommseq is ',ppm_ncommseq(topoid),lidir
!          CALL ppm_write(ppm_rank,'ppm_util_commopt_cart',mesg,info)
!          DO idir=1,ppm_ncommseq(topoid)
!             WRITE(mesg,*) ' entry ',ppm_icommseq(idir,topoid)
!             CALL ppm_write(ppm_rank,'ppm_util_commopt_cart',mesg,info)
!          END DO
!       END IF

      !  done
      !-----------------------------------------------------
      !-------------------------------------------------------------------------
      !  Everybody gets the memory needed
      !-------------------------------------------------------------------------
!       WRITE(mesg,*) 'ppm_icommseq ',ppm_icommseq(1:9,topoid)
!       CALL ppm_write(ppm_rank,'commopt',mesg,info)
!       WRITE(mesg,*) 'ppm_icommseq ',ppm_icommseq(10:18,topoid)
!       CALL ppm_write(ppm_rank,'commopt',mesg,info)
!       WRITE(mesg,*) 'ppm_icommseq ',ppm_icommseq(19:27,topoid)
!       CALL ppm_write(ppm_rank,'commopt',mesg,info)


      !-------------------------------------------------------------------------
      !  Mark this topology as done
      !-------------------------------------------------------------------------
      topo%isoptimized = .TRUE.

      !-------------------------------------------------------------------------
      !  Deallocate temporary storage
      !-------------------------------------------------------------------------
      iopt = ppm_param_dealloc
      CALL ppm_alloc(ichecked,ldu,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_util_commopt_cart',     &
              &            'number of comm.rounds PPM_NCOMMSEQ',__LINE__,info)
         GOTO 9999
      ENDIF

#endif

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_util_commopt_cart',t0,info)
      RETURN

      CONTAINS
      SUBROUTINE check
          CALL ppm_check_topoid(topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_util_commopt_cart',  &
     &            'topoid out of range',__LINE__,info)
              GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check

    END SUBROUTINE ppm_util_commopt_cart
