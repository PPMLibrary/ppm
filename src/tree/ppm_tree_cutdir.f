      !-------------------------------------------------------------------------
      !  Subroutine   :                  ppm_tree_cutdir
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

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_tree_cutdir_s(xp,Npart,weights,min_box,max_box,   &
     &    cutbox,ncut,fixed,minboxsize,icut,info,pcost)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_tree_cutdir_d(xp,Npart,weights,min_box,max_box,   &
     &    cutbox,ncut,fixed,minboxsize,icut,info,pcost)
#endif
      !!! This routine finds the best cutting directions and
      !!! positions for the given cut directions.
      !!!
      !!! [NOTE]
      !!! The cost of a particle is counted as pcost (if
      !!! given) or 1. For meshes, the cost is 1 per mesh
      !!! point. For the geometry part, the cost of a
      !!! box is given by its volume. Use weights(3) to
      !!! adjust this if needed.
      !!!
      !!! [NOTE]
      !!! The best direction cut is given by min cost of
      !!! geometry&mesh and max cost of moment of inertia
      !!! (particles cost)
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_tree
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_write
      IMPLICIT NONE
#if   __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single
#elif __KIND == __DOUBLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_double
#endif
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: xp
      !!! Particle positions
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: min_box
      !!! Minimum coordinate of the  boxes
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: max_box
      !!! Maximum coordinate of the  boxes
      REAL(MK), DIMENSION(3  ), INTENT(IN   ) :: minboxsize
      !!! Minimum size of a box in all directions.
      REAL(MK), DIMENSION(3  ), INTENT(IN   ) :: weights
      !!! Weights for the three cost contributions: particles, mesh,
      !!! geometry for finding the cut directions.
      REAL(MK), DIMENSION(:  ), OPTIONAL, INTENT(IN) :: pcost
      !!! Argument of length Npart, specifying the
      !!! computational cost of each particle.
      INTEGER                 , INTENT(IN   ) :: Npart
      !!! Number of particles
      INTEGER                 , INTENT(IN   ) :: ncut
      !!! Number of cut directions
      INTEGER                 , INTENT(IN   ) :: cutbox
      !!! ID of box to be cut
      LOGICAL , DIMENSION(:  ), INTENT(IN   ) :: fixed
      !!! Set to `TRUE` for dimensions which must not be cut
      INTEGER , DIMENSION(:  ), POINTER       :: icut
      !!! Directions of best cut. icut=i means: cutting plane is orthogonal
      !!! to i-th coordinate axis. index: 1..ncut. The directions are sorted
      !!! with the most favorable first.
      INTEGER                 , INTENT(  OUT) :: info
      !!! Return status, 0 on success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(ppm_dim)            :: len_box
      REAL(MK), DIMENSION(ppm_dim)            :: cgeom,cmesh,cpart,ctotal
      REAL(MK), DIMENSION(ppm_dim)            :: cpart_tot,shift
      LOGICAL , DIMENSION(ppm_dim)            :: cutable
      INTEGER , DIMENSION(ppm_dim)            :: isort
      REAL(MK), DIMENSION(ncut+1)             :: pc,pcsum
      REAL(MK)                                :: t0,dm
      REAL(MK)                                :: csum,csuminv,lmyeps
      REAL(MK)                                :: x,y,z,x2,y2,z2
      INTEGER                                 :: i,j,k,ip,cutdir,ncp1
#ifdef __MPI
      INTEGER                                 :: MPTYPE
#endif

      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_tree_cutdir',t0,info)
#if   __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif

      !-------------------------------------------------------------------------
      !  Check input arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug .GT. 0) THEN
        CALL check
        IF (info .NE. 0) GOTO 9999
      ENDIF
      
      !-------------------------------------------------------------------------
      !  If we have less than 1 direction to cut, we are done
      !-------------------------------------------------------------------------
      IF (ncut .LT. 1) THEN
          IF (ppm_debug .GT. 0) THEN
              CALL ppm_write(ppm_rank,'ppm_tree_cutdir',   &
     &            'No cut directions present. Exiting.',info)
          ENDIF
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Compute the extension of the box
      !-------------------------------------------------------------------------
      IF (ppm_dim .GT. 2) THEN
          len_box(1) = max_box(1,cutbox)-min_box(1,cutbox)
          len_box(2) = max_box(2,cutbox)-min_box(2,cutbox)
          len_box(3) = max_box(3,cutbox)-min_box(3,cutbox)
      ELSE
          len_box(1) = max_box(1,cutbox)-min_box(1,cutbox)
          len_box(2) = max_box(2,cutbox)-min_box(2,cutbox)
      ENDIF

      !-------------------------------------------------------------------------
      !  Determine cuttable directions of this box
      !-------------------------------------------------------------------------
      cutable = .FALSE.
      ip = 0
      DO i=1,ppm_dim
          IF ((.NOT.fixed(i)) .AND. (len_box(i)-(2.0_MK*minboxsize(i))   &
     &        .GT.lmyeps*len_box(i))) THEN
              ip = ip + 1
              cutable(i) = .TRUE.      
          ENDIF
      ENDDO

      !-------------------------------------------------------------------------
      !  Exit if box has not enough cuttable directions
      !-------------------------------------------------------------------------
      IF (ip .LT. ncut) THEN
          IF (ppm_debug .GT. 0) THEN
              CALL ppm_write(ppm_rank,'ppm_tree_cutdir',   &
     &            'Not enough cutable directions! Exiting.',info)
          ENDIF
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  If this is an octtree in 3d or a quad-tree in 2d it is trivial
      !-------------------------------------------------------------------------
      IF (ncut .EQ. 3 .AND. ppm_dim .EQ. 3) THEN
          icut(1) = 1
          icut(2) = 2
          icut(3) = 3
          GOTO 9999
      ELSEIF (ncut .EQ. 2 .AND. ppm_dim .EQ. 2) THEN
          icut(1) = 1
          icut(2) = 2
          GOTO 9999
      ENDIF

           
#ifdef __MPI
      !-------------------------------------------------------------------------
      !  Determine MPI data type
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      MPTYPE = MPI_REAL
#elif __KIND == __DOUBLE_PRECISION 
      MPTYPE = MPI_DOUBLE_PRECISION
#endif
#endif
       
      !------------------------------------------------------------------------
      !  Compute the geometry cost:
      !  - for 3D the cgeom(1) = len_box(2)*len_box(3)
      !               cgeom(2) = len_box(1)*len_box(3)
      !               cgeom(3) = len_box(1)*len_box(2)
      !  In order to obtain the sum of the geometry costs on all axes 1, 
      !  the cost on each axis is normalized with the sum of all costs.  
      !------------------------------------------------------------------------
      ctotal = 0.0_MK
      cgeom  = 0.0_MK
      IF (weights(3) .NE. 0.0_MK) THEN
          IF (ppm_dim .GT. 2) THEN
              cgeom(1) = len_box(2)*len_box(3)
              cgeom(2) = len_box(1)*len_box(3)
              cgeom(3) = len_box(1)*len_box(2)
              ! normalization: total cost contribution of each part needs to
              ! sum up to 1 in order to get correct weighting.
              csum     = cgeom(1)+cgeom(2)+cgeom(3)
              csuminv  = 1.0_MK/csum
              cgeom(1) = cgeom(1)*csuminv
              cgeom(2) = cgeom(2)*csuminv
              cgeom(3) = cgeom(3)*csuminv
              ctotal(1)= weights(3)*cgeom(1)
              ctotal(2)= weights(3)*cgeom(2)
              ctotal(3)= weights(3)*cgeom(3)
          ELSE
              cgeom(1) = len_box(2)
              cgeom(2) = len_box(1)
              ! normalization: total cost contribution of each part needs to
              ! sum up to 1 in order to get correct weighting.
              csum     = cgeom(1)+cgeom(2)
              csuminv  = 1.0_MK/csum
              cgeom(1) = cgeom(1)*csuminv
              cgeom(2) = cgeom(2)*csuminv
              ctotal(1)= weights(3)*cgeom(1)
              ctotal(2)= weights(3)*cgeom(2)
          ENDIF
      ENDIF

      !------------------------------------------------------------------------
      !  Compute the mesh cost: 
      !  - for 2D the cmesh(i) Nm_box(i,cutbox)
      !  - for 3D the cmesh(1) = Nm_box(2,cutbox)*Nm_box(3,cutbox)
      !               cmesh(2) = Nm_box(1,cutbox)*Nm_box(3,cutbox)
      !               cmesh(3) = Nm_box(1,cutbox)*Nm_box(2,cutbox)
      !------------------------------------------------------------------------
      cmesh = 0.0_MK
      IF (have_mesh .AND. weights(2) .NE. 0.0_MK) THEN
          IF (ppm_dim .GT. 2) THEN
              cmesh(1) = Nm_box(2,cutbox)*Nm_box(3,cutbox)
              cmesh(2) = Nm_box(1,cutbox)*Nm_box(3,cutbox)
              cmesh(3) = Nm_box(1,cutbox)*Nm_box(2,cutbox) 
              ! normalization: total cost contribution of each part needs to
              ! sum up to 1 in order to get correct weighting.
              csum     = REAL(cmesh(1),MK)+REAL(cmesh(2),MK)+REAL(cmesh(3),MK)
              csuminv  = 1.0_MK/csum
              cmesh(1) = cmesh(1)*csuminv
              cmesh(2) = cmesh(2)*csuminv
              cmesh(3) = cmesh(3)*csuminv
              ctotal(1)= ctotal(1) + (weights(2)*cmesh(1))
              ctotal(2)= ctotal(2) + (weights(2)*cmesh(2))
              ctotal(3)= ctotal(3) + (weights(2)*cmesh(3))
          ELSE
              cmesh(1) = Nm_box(2,cutbox)
              cmesh(2) = Nm_box(1,cutbox)
              ! normalization: total cost contribution of each part needs to
              ! sum up to 1 in order to get correct weighting.
              csum     = REAL(cmesh(1),MK)+REAL(cmesh(2),MK)
              csuminv  = 1.0_MK/csum
              cmesh(1) = cmesh(1)*csuminv
              cmesh(2) = cmesh(2)*csuminv
              ctotal(1)= ctotal(1) + (weights(2)*cmesh(1))
              ctotal(2)= ctotal(2) + (weights(2)*cmesh(2))
          ENDIF
      ENDIF
      
      !------------------------------------------------------------------------
      !  Compute the particles cost: the moments of inertia
      !------------------------------------------------------------------------
      cpart = 0.0_MK
      IF (have_particles .AND. weights(1) .NE. 0.0_MK) THEN      
          !---------------------------------------------------------------------
          !  Compute the cost of local particles
          !---------------------------------------------------------------------
          shift(1) = 0.5_MK*(min_box(1,cutbox)+max_box(1,cutbox))
          shift(2) = 0.5_MK*(min_box(2,cutbox)+max_box(2,cutbox))
          IF (ppm_dim .EQ. 2) THEN
              DO k=tree_lhbx(1,cutbox),tree_lhbx(2,cutbox)
                  ip = tree_lpdx(k)
                  x = xp(1,ip)-shift(1)
                  y = xp(2,ip)-shift(2)
                  x2 = x*x
                  y2 = y*y
                  dm = 1.0_MK
                  IF (PRESENT(pcost)) dm = pcost(ip)
                  cpart(1) = cpart(1) + (y2*dm)
                  cpart(2) = cpart(2) + (x2*dm)
              ENDDO
          ELSE
              shift(3) = 0.5_MK*(min_box(3,cutbox)+max_box(3,cutbox))  
              DO k=tree_lhbx(1,cutbox),tree_lhbx(2,cutbox)
                  ip = tree_lpdx(k)
                  x = xp(1,ip)-shift(1)
                  y = xp(2,ip)-shift(2)
                  z = xp(3,ip)-shift(3)
                  x2 = y*y+z*z
                  y2 = x*x+z*z
                  z2 = x*x+y*y
                  dm = 1.0_MK
                  IF (PRESENT(pcost)) dm = pcost(ip)
                  cpart(1) = cpart(1) + (x2*dm) 
                  cpart(2) = cpart(2) + (y2*dm)
                  cpart(3) = cpart(3) + (z2*dm)
              ENDDO
          ENDIF
#ifdef __MPI
          !---------------------------------------------------------------------
          !  Allreduce of particle costs 
          !---------------------------------------------------------------------
          CALL MPI_Allreduce(cpart,cpart_tot,ppm_dim,MPTYPE,MPI_SUM,    &
     &        ppm_comm,info)
          IF (info.NE.0) THEN
              info = ppm_error_fatal
              CALL ppm_error(ppm_err_mpi_fail,'ppm_tree_cutdir',   &
     &            'MPI_Allreduce of inertia tensor failed!',__LINE__,info)
              GOTO 9999
          ENDIF 
          cpart = cpart_tot
#endif
         
          !---------------------------------------------------------------------
          !  Normalize the particle cost
          !---------------------------------------------------------------------
          IF (ppm_dim .GT. 2) THEN
              csum      = cpart(1)+cpart(2)+cpart(3)
              csuminv   = 1.0_MK/csum
              cpart(1)  = cpart(1)*csuminv
              cpart(2)  = cpart(2)*csuminv
              cpart(3)  = cpart(3)*csuminv
              ctotal(1) = ctotal(1)+(weights(1)*cpart(1))           
              ctotal(2) = ctotal(2)+(weights(1)*cpart(2))           
              ctotal(3) = ctotal(3)+(weights(1)*cpart(3))           
          ELSE
              csum      = cpart(1)+cpart(2)
              csuminv   = 1.0_MK/csum
              cpart(1)  = cpart(1)*csuminv
              cpart(2)  = cpart(2)*csuminv
              ctotal(1) = ctotal(1)+(weights(1)*cpart(1))           
              ctotal(2) = ctotal(2)+(weights(1)*cpart(2))           
          ENDIF
      ENDIF ! have_particles
   
      !-------------------------------------------------------------------------
      !  Sort directions in decreasing order by total cost. We do not need to 
      !  normalize ctotal by SUM(weights), since we are only interested in the 
      !  ORDER of directions and not the actual cost values.
      !-------------------------------------------------------------------------
      IF (ppm_dim .EQ. 2) THEN
          isort(1:2) = (/1,2/)
          IF (ctotal(2) .LT. ctotal(1)) isort(1:2) = (/2,1/)
      ELSE
          isort(1:3) = (/1,2,3/)
          IF (ctotal(2) .LT. ctotal(1)) isort(1:3) = (/2,1,3/)
          IF (ctotal(3).LT.ctotal(1).AND.ctotal(3).LT.ctotal(2)) THEN 
              isort(3) = isort(2)
              isort(2) = isort(1)
              isort(1) = 3
          ELSEIF (ctotal(3).LT.ctotal(isort(2))) THEN
              isort(3) = isort(2)
              isort(2) = 3
          ENDIF
      ENDIF

      !-------------------------------------------------------------------------
      !  Assign ncut best non-fixed directions in ascending order to icut
      !-------------------------------------------------------------------------
      ip = 0
      i  = 1
      DO WHILE ((ip .LT. ncut) .AND. (i .LE. ppm_dim))
          k = isort(i)
          IF (cutable(k)) THEN
              ip = ip + 1
              icut(ip) = k
          ENDIF
          i = i + 1
      ENDDO

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_tree_cutdir',t0,info)
      RETURN
      CONTAINS
      SUBROUTINE check
         IF (cutbox .LE. 0) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_tree_cutdir',     &
     &          'cutbox must be > 0 !',__LINE__,info)
            GOTO 8888
         ENDIF
         IF (SIZE(min_box,2) .LT. cutbox) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_tree_cutdir',     &
     &          'size of min_box must be at least cutbox !',__LINE__,info)
            GOTO 8888
         ENDIF
         IF (SIZE(max_box,2) .LT. cutbox) THEN
            info = ppm_error_error
            CALL ppm_error(ppm_err_argument,'ppm_tree_cutdir',     &
     &          'size of max_box must be at least cutbox !',__LINE__,info)
            GOTO 8888
         ENDIF
         DO i=1,ppm_dim
            IF (min_box(i,cutbox) .GT. max_box(i,cutbox)) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_tree_cutdir',     &
     &             'min_box must be <= max_box !',__LINE__,info)
               GOTO 8888
            ENDIF
         ENDDO
 8888    CONTINUE
      END SUBROUTINE check
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE ppm_tree_cutdir_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_tree_cutdir_d
#endif
 
