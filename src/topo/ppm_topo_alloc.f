      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_topo_alloc
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

      SUBROUTINE ppm_topo_alloc(topoid,nsubs,nsublist,maxneigh,prec,info)
      !!! This routine (re-)allocates a topology object and all its members

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_topo_typedef
      IMPLICIT NONE
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER, INTENT(INOUT) :: topoid
      !!! Topology id structure to be (re)allocated
      !!!
      !!! if topoid == 0 then a new topology should be allocated, else the
      !!! topology with id == topoid is reallocated
      INTEGER, INTENT(IN   ) :: nsubs
      !!! Total number of subs on all procs
      INTEGER, INTENT(IN   ) :: nsublist
      !!! Local number of subs on this proc
      INTEGER, INTENT(IN   ) :: maxneigh
      !!! Maximum number of neighbours of any sub on this processor
      INTEGER, INTENT(IN   ) :: prec
      !!! Precision for storage. One of:
      !!!
      !!! * ppm_kind_single
      !!! * ppm_kind_double
      INTEGER, INTENT(  OUT) :: info
      !!! Return status, on success 0.


      !-------------------------------------------------------------------------
      !  Local variables
      !-------------------------------------------------------------------------
      TYPE(ppm_t_topo), POINTER :: topo

      TYPE(ppm_t_ptr_topo), DIMENSION(:), POINTER :: temptopo

      REAL(ppm_kind_double) :: t0

      INTEGER, DIMENSION(3) :: ldc
      !!! Number of elements in all dimensions for allocation
      INTEGER               :: iopt
      !!! allocation mode, see one of ppm_alloc_* subroutines.
      INTEGER               :: nsubmax
      !!! contains the MAX of topo%nsubs and 1
      INTEGER               :: nsublistmax
      !!! contains the MAX of topo%nsublist and 1
      INTEGER               :: topo_idx
      !!! local variable holding the index in the ppm_topo array with the
      !!! topology that should be (re)allocated
      INTEGER               :: i

      CHARACTER(LEN=ppm_char) :: caller="ppm_topo_alloc"

      !-------------------------------------------------------------------------
      !  Initialize
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN
         CALL check
         IF (info.NE.0) GOTO 9999
      ENDIF

      IF (topoid.EQ.0) THEN
         !-----------------------------------------------------------------------
         ! We have to create a new topology
         !-----------------------------------------------------------------------
         IF (ppm_next_avail_topo.EQ.ppm_param_undefined) THEN
            ! This is the first topology, initialize everything
            ALLOCATE(ppm_topo(8),STAT=info)
            or_fail_alloc("Could not allocate ppm_topo",ppm_error=ppm_error_fatal)

            ! make sure all pointers are nullified
            DO i=1,8
               NULLIFY(ppm_topo(i)%t)
            ENDDO

            ppm_next_avail_topo = 1

            topoid = ppm_next_avail_topo

            ppm_next_avail_topo = ppm_next_avail_topo + 1
         ELSE IF (ppm_next_avail_topo.GT.SIZE(ppm_topo)) THEN
            ! We need more space in the ppm_topo array, enlarge
            ALLOCATE(temptopo(SIZE(ppm_topo)),STAT=info)
            or_fail_alloc("Could not allocate temptopo")

            DO i=1,SIZE(temptopo)
               temptopo(i) = ppm_topo(i)
            ENDDO

            DEALLOCATE(ppm_topo,STAT=info)
            or_fail_dealloc("Could not deallocate ppm_topo")

            ALLOCATE(ppm_topo(2*SIZE(temptopo)),STAT=info)
            or_fail_alloc("Could not reallocate ppm_topo")

            DO i=1,SIZE(temptopo)
               ppm_topo(i) = temptopo(i)
            ENDDO

            DO i=SIZE(temptopo)+1,SIZE(ppm_topo)
               NULLIFY(ppm_topo(i)%t)
            ENDDO

            DEALLOCATE(temptopo,STAT=info)
            or_fail_dealloc("Could not deallocate temptopo")
            NULLIFY(temptopo)

            topoid = ppm_next_avail_topo

            ppm_next_avail_topo = ppm_next_avail_topo + 1
         ELSE
            ! We are using this topoid for the new topology
            topoid = ppm_next_avail_topo

            ! We need to find the next available topoid
            ppm_next_avail_topo=-1

            DO i=1,SIZE(ppm_topo)
               IF (i.EQ.topoid) CYCLE

               IF (ASSOCIATED(ppm_topo(i)%t)) THEN
                  IF (ASSOCIATED(ppm_topo(i)%t%bcdef)) THEN
                     CYCLE
                  ELSE
                     ppm_next_avail_topo = i
                     EXIT
                  ENDIF
               ELSE
                  ppm_next_avail_topo = i
                  EXIT
               ENDIF
            ENDDO

            ppm_next_avail_topo=MERGE(ppm_next_avail_topo,SIZE(ppm_topo)+1,ppm_next_avail_topo.GT.0)
         ENDIF
      ENDIF

      CALL ppm_alloc(ppm_topo(topoid)%t,ppm_param_alloc_fit,info)
      or_fail_alloc("Could not allocate ppm_topo elements")

      topo => ppm_topo(topoid)%t

      topo%ID = topoid

      !-------------------------------------------------------------------------
      !  Allocate memory for the domain
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = ppm_dim
      SELECT CASE (prec)
      CASE (ppm_kind_single)
         CALL ppm_alloc(topo%min_physs,ldc,iopt,info)
         or_fail_alloc("min extent of domain TOPO%MIN_PHYSS",ppm_error=ppm_error_fatal)

         CALL ppm_alloc(topo%max_physs,ldc,iopt,info)
         or_fail_alloc("max extent of domain TOPO%MAX_PHYSS",ppm_error=ppm_error_fatal)
      CASE (ppm_kind_double)
          CALL ppm_alloc(topo%min_physd,ldc,iopt,info)
          or_fail_alloc("min extent of domain TOPO%MIN_PHYSD",ppm_error=ppm_error_fatal)

          CALL ppm_alloc(topo%max_physd,ldc,iopt,info)
          or_fail_alloc("max extent of domain TOPO%MAX_PHYSD",ppm_error=ppm_error_fatal)
      END SELECT
      !-------------------------------------------------------------------------
      !  The MAX of naublist and 1 is needed to avoid allocation failures
      !  if a processor has 0 subs, same goes for nsubs
      !-------------------------------------------------------------------------
      nsubmax = MAX(nsubs,1)
      nsublistmax = MAX(nsublist,1)

      !-------------------------------------------------------------------------
      !  Allocate memory for the subdomains
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = ppm_dim
      ldc(2) = nsubs
      SELECT CASE (prec)
      CASE (ppm_kind_single)
          CALL ppm_alloc(topo%min_subs,ldc,iopt,info)
          or_fail_alloc("min extent of subs TOPO%MIN_SUBS",ppm_error=ppm_error_fatal)

          CALL ppm_alloc(topo%max_subs,ldc,iopt,info)
          or_fail_alloc("max extent of subs TOPO%MAX_SUBS",ppm_error=ppm_error_fatal)
      CASE (ppm_kind_double)
          CALL ppm_alloc(topo%min_subd,ldc,iopt,info)
          or_fail_alloc("min extent of subs TOPO%MIN_SUBD",ppm_error=ppm_error_fatal)

          CALL ppm_alloc(topo%max_subd,ldc,iopt,info)
          or_fail_alloc("max extent of subs TOPO%MAX_SUBD",ppm_error=ppm_error_fatal)
      END SELECT

      !-------------------------------------------------------------------------
      !  Allocate memory for the sub-to-proc mapping
      !-------------------------------------------------------------------------
      ldc(1) = nsubs
      CALL ppm_alloc(topo%sub2proc,ldc,iopt,info)
      or_fail_alloc("global subs to proc map TOPO%SUB2PROC",ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Allocate memory for list of local subs
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = nsubmax
      CALL ppm_alloc(topo%isublist,ldc,iopt,info)
      or_fail_alloc("local sub list TOPO%ISUBLIST",ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Allocate memory for the boundary conditions
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = 2*ppm_dim
      CALL ppm_alloc(topo%bcdef,ldc,iopt,info)
      or_fail_alloc("boundary conditions TOPO%BCDEF",ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Allocate memory for the subdomain costs
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = nsubs
      SELECT CASE (prec)
      CASE (ppm_kind_single)
         CALL ppm_alloc(topo%sub_costs,ldc,iopt,info)
      CASE (ppm_kind_double)
         CALL ppm_alloc(topo%sub_costd,ldc,iopt,info)
      END SELECT
      or_fail_alloc("subdomain costs TOPO%SUB_COST",ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Allocate memory for the numbers of neighbors of each local sub on
      !  the current processor i.e, a subset of the nneigh(1:nsubs) list.
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = nsublistmax
      CALL ppm_alloc(topo%nneighsubs,ldc,iopt,info)
      or_fail_alloc("number of neighbors of local subs TOPO%NNEIGHSUBS",ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Allocate memory for the global IDs of all the neighbors.
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = maxneigh
      ldc(2) = nsublist
      CALL ppm_alloc(topo%ineighsubs,ldc,iopt,info)
      or_fail_alloc("neighbors of local subs TOPO%INEIGHSUBS",ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Allocate memory for BC for the subs on the current processor
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = 2*ppm_dim
      ldc(2) = nsubmax
      CALL ppm_alloc(topo%subs_bc,ldc,iopt,info)
      or_fail_alloc("BCs for subs on local processor TOPO%SUBS_BC",ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Allocate memory for the list of neighboring processors.
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = 26    ! just guessing :-) FIXME: Why?
      CALL ppm_alloc(topo%ineighproc,ldc,iopt,info)
      or_fail_alloc("list of neighbor processors TOPO%INEIGHPROC",ppm_error=ppm_error_fatal)

      !-------------------------------------------------------------------------
      !  Initialize optmization status flags
      !-------------------------------------------------------------------------
      topo%isoptimized = .FALSE.

      !-------------------------------------------------------------------------
      !  Mark this topology as not defined
      !-------------------------------------------------------------------------
      topo%isdefined = .FALSE.

      !-------------------------------------------------------------------------
      !  Set the precision for this topology
      !-------------------------------------------------------------------------
      topo%prec = prec

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
      9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
          IF (nsubs.LE.0) THEN
             fail("nsubs must be >0",exit_point=8888)
          ENDIF
          IF (nsublist.LT.0) THEN
             fail("nsublist is not >=0",exit_point=8888,ppm_error=ppm_error_notice)
          ENDIF
          IF (nsubs.LT.nsublist) THEN
             fail("total number of subs is smaller than local",exit_point=8888)
          ENDIF
          IF (maxneigh.LT.0) THEN
             fail("maxneigh must be >=0",exit_point=8888)
          ENDIF
          IF (topoid.NE.0) THEN
             IF ((topoid.GT.SIZE(ppm_topo)).OR.(topoid.LT.1)) THEN
                fail("topoid indexing outside ppm_topo",exit_point=8888)
              ENDIF
          ENDIF
      8888 CONTINUE
      END SUBROUTINE check
      END SUBROUTINE ppm_topo_alloc
