      !-------------------------------------------------------------------------
      !  Subroutine   :                   ppm_topo_get_meshinfo
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

      SUBROUTINE ppm_topo_get_meshinfo(topoid,meshid,nm,istart,ndata,maxndata,&
      &                                isublist,nsublist,info)
      !!! This routine returns the subdomain boundaries and boundary
      !!! conditions

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_typedef
      USE ppm_module_alloc
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER                         , INTENT(IN   ) :: topoid
      !!! Topology ID
      INTEGER                         , INTENT(IN   ) :: meshid
      !!! Mesh ID
      INTEGER, DIMENSION(:  ), POINTER                :: nm
      !!! global number of mesh points in computational domain
      INTEGER, DIMENSION(:,:), POINTER                :: istart
      !!! Start indices (i,j[,k]) (first index) of sub mesh isub (second
      !!! index) in *global* mesh.
      !!!
      !!! [CAUTION]
      !!! This array holds the starting indeces for *all*
      !!! subdomains in the *whole* physical domain (i.e. all processors).
      !!! Use `topo%isublist` to find the indeces for the grids on the
      !!! processor
      INTEGER, DIMENSION(:,:), POINTER                :: ndata
      !!! Returns number of grid POINTS in x,y[,z] (first index) of sub mesh
      !!! isub (second index). Includes the points ON the sub boundaries.
      !!!
      !!! [CAUTION]
      !!! This array holds the number of points for *all* subdomains in the
      !!! *whole* physical domain (i.e. all processors).
      !!! Use `topo%isublist` to find the indeces for the grids on the
      !!! processor
      INTEGER, DIMENSION(ppm_dim)      , INTENT(  OUT) :: maxndata
      !!! Returns maximum number of grid points over all subdomains
      INTEGER, DIMENSION(:  ), POINTER                 :: isublist
      !!! list of subs of the current processor.
      !!! 1st index: local sub number.
      INTEGER                         ,  INTENT(  OUT) :: nsublist
      !!! number of subs on the current processor.
      INTEGER                          , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)          :: t0
      INTEGER                        :: i,isub,d
      LOGICAL                        :: valid
      TYPE(ppm_t_topo),      POINTER :: topo
      TYPE(ppm_t_equi_mesh), POINTER :: mesh
      INTEGER, DIMENSION(2)          :: lda
      INTEGER                        :: iopt
      !-------------------------------------------------------------------------
      !  Externals
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart('ppm_topo_get_meshinfo',t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF


      topo => ppm_topo(topoid)%t
      mesh => topo%mesh(meshid)

      !-------------------------------------------------------------------------
      !  Allocate memory for structures
      !-------------------------------------------------------------------------

      iopt = ppm_param_alloc_fit
      lda(1) = ppm_dim
      lda(2) = topo%nsubs
      CALL ppm_alloc(nm,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_get_meshinfo',    &
     &        'global mesh size NM',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(istart,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_get_meshinfo',    &
     &        'sub mesh start indices ISTART',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ndata,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_get_meshinfo',    &
     &        'sub mesh sizes NDATA',__LINE__,info)
          GOTO 9999
      ENDIF

      lda(1) = topo%nsublist
      CALL ppm_alloc(isublist,lda,iopt,info)
      IF (info .NE. 0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_topo_get_meshinfo',    &
     &        'sub mesh sizes ISUBLIST',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Copy data
      !-------------------------------------------------------------------------

      nsublist = topo%nsublist
      FORALL (i=1:nsublist) isublist(i) = topo%isublist(i)
      nm = mesh%Nm
      FORALL(i=1:topo%nsubs, d=1:ppm_dim)
          ndata(d,i) = mesh%nnodes(d,i)
          istart(d,i) = mesh%istart(d,i)
      END FORALL
      maxndata = 0
      DO i=1,nsublist
          isub = isublist(i)
          DO d=1,ppm_dim
              IF (ndata(d,isub).GT.maxndata(d)) maxndata(d) = ndata(d,isub)
          ENDDO
      ENDDO

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_topo_get_meshinfo',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF (.NOT. ppm_initialized) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_ppm_noinit,'ppm_topo_get',  &
     &            'Please call ppm_init first!',__LINE__,info)
              GOTO 8888
          ENDIF
          CALL ppm_check_topoid(topoid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_get_meshinfo',  &
     &            'Topology ID is invalid!',__LINE__,info)
              GOTO 8888
          ENDIF
          CALL ppm_check_meshid(topoid,meshid,valid,info)
          IF (.NOT. valid) THEN
              info = ppm_error_error
              CALL ppm_error(ppm_err_argument,'ppm_topo_get_meshinfo',  &
     &            'Mesh ID is invalid!',__LINE__,info)
              GOTO 8888
          ENDIF
 8888     CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_topo_get_meshinfo
