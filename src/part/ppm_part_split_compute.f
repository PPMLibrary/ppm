      !-------------------------------------------------------------------------
      !  Subroutine   :               ppm_part_split_compute
      !-------------------------------------------------------------------------
      ! Copyright (c) 2010 CSE Lab (ETH Zurich), MOSAIC Group (ETH Zurich), 
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

#if    __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_part_split_compute_s(topoid,xp,Nall,isymm,ghostsize,info)
#elif  __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_part_split_compute_d(topoid,xp,Nall,isymm,ghostsize,info)
#endif 
      !!! This routine splits a set of new particles into 
      !!! real particles and ghost particles.
      !!! ==============================================================


      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_typedef
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
      USE ppm_module_write
      USE ppm_module_check_id
      USE ppm_module_util_commopt
      USE ppm_module_data_buffers_add
      IMPLICIT NONE
#if    __KIND == __SINGLE_PRECISION  | __KIND_AUX == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#else
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(IN   ) :: topoid
      !!! ID of current topology
      REAL(MK), DIMENSION(:,:), INTENT(IN   ),POINTER :: xp
      !!! The position of the new particles to be added
      INTEGER                 , INTENT(IN   ) :: Nall
      !!! The number of new particles (on the local processor)
      INTEGER                 , INTENT(IN   ) :: isymm
      !!! Indicator for the use of symmetry
      !!!
      !!! * isymm > 0 use symmetry
      !!! * isymm = 0 do not use symmetry
      REAL(MK)                , INTENT(IN   ) :: ghostsize
      !!! The size of the ghost layer
      INTEGER                 , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER, DIMENSION(3)         :: ldu
      INTEGER                       :: i,j,k,isub
      INTEGER                       :: nlist1,nlist2,ipart,iopt
      REAL(MK)                      :: xminf,yminf,zminf ! full domain
      REAL(MK)                      :: xmaxf,ymaxf,zmaxf ! full domain
      REAL(MK)                      :: xmini,ymini,zmini ! inner domain
      REAL(MK)                      :: xmaxi,ymaxi,zmaxi ! inner domain
      REAL(MK), DIMENSION(ppm_dim)  :: min_phys
      REAL(MK), DIMENSION(ppm_dim)  :: max_phys
      REAL(MK), DIMENSION(ppm_dim)  :: len_phys
      REAL(MK)                      :: t0
      TYPE(ppm_t_topo),POINTER      :: topo => NULL()
      CHARACTER(LEN=ppm_char)       :: caller = 'ppm_part_split_compute'
      INTEGER, DIMENSION(2*ppm_dim) :: ibc

      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart(caller,t0,info)

      topo => ppm_topo(topoid)%t

      !-------------------------------------------------------------------------
      !  Compute the size of the computational box on this topology
      !-------------------------------------------------------------------------
#if __KIND == __DOUBLE_PRECISION
      min_phys(:) = topo%min_physd(:)
      max_phys(:) = topo%max_physd(:)
#else
      min_phys(:) = topo%min_physs(:)
      max_phys(:) = topo%max_physs(:)
#endif
      len_phys(:) = max_phys(:) - min_phys(:)

      !-------------------------------------------------------------------------
      !  Allocate memory for the list of particle on the local processor that
      !  may be ghosts on other processors
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_grow
      ldu(1) = Nall
      CALL ppm_alloc(ilist1,ldu,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,caller,'list1',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ilist2,ldu,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,caller,'list2',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(ighost,ldu,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,caller,'ighost',__LINE__,info)
          GOTO 9999
      ENDIF

      ldu(1) = Nall
      CALL ppm_alloc(modify%idx_real_new,ldu,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,caller,'XPR',__LINE__,info)
          GOTO 9999
      ENDIF
      CALL ppm_alloc(modify%idx_ghost_new,ldu,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,caller,'XPG',__LINE__,info)
          GOTO 9999
      ENDIF
      modify%Nrnew = 0
      modify%Ngnew = 0

      !-------------------------------------------------------------------------
      !  List ilist1() holds the particles we are currently considering 
      !  List ilist2() holds the particles that have not yet been associated
      !  with a sub
      !-------------------------------------------------------------------------
      DO i=1,Nall
         ilist1(i) = i
      ENDDO
      nlist1  = Nall
      nlist2  = 0
       
      ibc(:) = 0
      DO i=1,2*ppm_dim
          SELECT CASE (topo%bcdef(i))
          CASE (ppm_param_bcdef_periodic)
              ibc(i) = ppm_param_bcdef_periodic
          CASE (ppm_param_bcdef_symmetry)
              ibc(i) = ppm_param_bcdef_symmetry
          CASE (ppm_param_bcdef_antisymmetry)
              ibc(i) = ppm_param_bcdef_symmetry
          CASE (ppm_param_bcdef_neumann)
              ibc(i) = ppm_param_bcdef_symmetry
          CASE (ppm_param_bcdef_dirichlet)
              ibc(i) = ppm_param_bcdef_symmetry
          END SELECT
      ENDDO

      IF (ibc(1).EQ.ppm_param_bcdef_periodic) THEN
#include "../map/ghost_map_periodic_bc_x_add.inc"
      ELSEIF (ibc(1).EQ.ppm_param_bcdef_symmetry) THEN
#include "../map/ghost_map_symmetry_bc_lx_add.inc"
      ENDIF  
      
      IF (ibc(2).EQ.ppm_param_bcdef_symmetry) THEN
#include "../map/ghost_map_symmetry_bc_ux_add.inc"
      ENDIF 
      !----------------------------------------------------------------------
      !  handle boundary conditions in y
      !----------------------------------------------------------------------
      IF (ibc(3).EQ.ppm_param_bcdef_periodic) THEN
#include "../map/ghost_map_periodic_bc_y_add.inc"
      ELSEIF (ibc(3).EQ.ppm_param_bcdef_symmetry) THEN
#include "../map/ghost_map_symmetry_bc_ly_add.inc"
      ENDIF
      IF (ibc(4).EQ.ppm_param_bcdef_symmetry) THEN
#include "../map/ghost_map_symmetry_bc_uy_add.inc"
      ENDIF
      !----------------------------------------------------------------------
      !  handle periodicity in z (if 3D)
      !----------------------------------------------------------------------
      IF (ppm_dim.EQ.3) THEN
         !-------------------------------------------------------------------
         !  yes, we split the if in two, since we do not know in what order
         !  the compiler will check and ppm_bcdef will only be allocated to
         !  four (4) in 2D
         !-------------------------------------------------------------------
         IF (ibc(5).EQ.ppm_param_bcdef_periodic) THEN
#include "../map/ghost_map_periodic_bc_z_add.inc"
         ELSEIF (ibc(5).EQ.ppm_param_bcdef_symmetry) THEN
#include "../map/ghost_map_symmetry_bc_lz_add.inc"
         ENDIF
         IF (ibc(6).EQ.ppm_param_bcdef_symmetry) THEN
#include "../map/ghost_map_symmetry_bc_uz_add.inc"
         ENDIF

      ENDIF ! of 3D


      !-------------------------------------------------------------------------
      !  Loop over the subs belonging to this processor and check if 
      !  the particles are within the sub. If this is the case, the particle is 
      !  added to ilist2.
      !  Once all subs have been considered, the remaining particles are ghosts.
      !-------------------------------------------------------------------------
      DO k=1,topo%nsublist
         !----------------------------------------------------------------------
         !  Initialize the second list counter to zero
         !----------------------------------------------------------------------
         nlist2 = 0

         !----------------------------------------------------------------------
         !  Get the (global) id of the sub
         !----------------------------------------------------------------------
         isub = topo%isublist(k)

         !----------------------------------------------------------------------
         !  Store the full extend of the sub in ?minf and ?maxf
         !----------------------------------------------------------------------
#if __KIND == __DOUBLE_PRECISION
         xminf = topo%min_subd(1,isub)
         xmaxf = topo%max_subd(1,isub)

         yminf = topo%min_subd(2,isub)
         ymaxf = topo%max_subd(2,isub)

         IF (ppm_dim.EQ.3) THEN
            zminf = topo%min_subd(3,isub)
            zmaxf = topo%max_subd(3,isub)
         ENDIF 
#else
         xminf = topo%min_subs(1,isub)
         xmaxf = topo%max_subs(1,isub)

         yminf = topo%min_subs(2,isub)
         ymaxf = topo%max_subs(2,isub)
         IF (ppm_dim.EQ.3) THEN
            zminf = topo%min_subs(3,isub)
            zmaxf = topo%max_subs(3,isub)
         ENDIF
#endif
         !----------------------------------------------------------------------
         !  loop over the remaining particles
         !----------------------------------------------------------------------
         IF (ppm_dim.EQ.2) THEN
            DO j=1,nlist1
               !----------------------------------------------------------------
               !  get the particle index
               !----------------------------------------------------------------
               ipart = ilist1(j)
   
               !----------------------------------------------------------------
               !  check if the particles belongs to this sub
               !----------------------------------------------------------------
               IF (xp(1,ipart).GE.xminf.AND.xp(1,ipart).LT.xmaxf.AND. &
     &             xp(2,ipart).GE.yminf.AND.xp(2,ipart).LT.ymaxf) THEN
                  !-------------------------------------------------------------
                  !  if yes, add this particle to the list of real particles
                  !-------------------------------------------------------------
                  modify%Nrnew     = modify%Nrnew + 1
                  modify%idx_real_new(modify%Nrnew) = ipart
               ELSE    
                  !-------------------------------------------------------------
                  !  If not on this sub we need to consider it further
                  !-------------------------------------------------------------
                  nlist2         = nlist2 + 1
                  ilist2(nlist2) = ipart
               ENDIF
            ENDDO
         ELSE
            DO j=1,nlist1
               !----------------------------------------------------------------
               !  get the particle index
               !----------------------------------------------------------------
               ipart = ilist1(j)
   
               !----------------------------------------------------------------
               !  check if the particles belongs to this sub
               !----------------------------------------------------------------
               IF (xp(1,ipart).GE.xminf.AND.xp(1,ipart).LT.xmaxf.AND. &
     &             xp(2,ipart).GE.yminf.AND.xp(2,ipart).LT.ymaxf.AND. &
     &             xp(3,ipart).GE.zminf.AND.xp(3,ipart).LT.zmaxf) THEN
                  !-------------------------------------------------------------
                  !  if yes, add this particle to the list of real particles
                  !-------------------------------------------------------------
                  modify%Nrnew     = modify%Nrnew + 1
                  modify%idx_real_new(modify%Nrnew) = ipart
               ELSE    
                  !-------------------------------------------------------------
                  !  If not on this sub we need to consider it further
                  !-------------------------------------------------------------
                  nlist2         = nlist2 + 1
                  ilist2(nlist2) = ipart
               ENDIF
            ENDDO
         ENDIF 

         !----------------------------------------------------------------------
         !  swap the lists
         !----------------------------------------------------------------------
         DO j=1,nlist2
            ilist1(j) = ilist2(j)
         ENDDO
         nlist1 = nlist2

      ENDDO ! end of subs on local processor


      !-------------------------------------------------------------------------
      !  At the end the nlist2 is the number of unassigned particles.
      !  These are ghosts.
      !-------------------------------------------------------------------------
      DO j=1,nlist2
          modify%idx_ghost_new(j) = ilist2(j)
      ENDDO
      modify%Ngnew                = nlist2


      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop(caller,t0,info)
      RETURN
#if    __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_part_split_compute_s
#elif  __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_part_split_compute_d
#endif
