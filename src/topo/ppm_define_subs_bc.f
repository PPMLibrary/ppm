      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_define_subs_bc
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
      SUBROUTINE define_subsbc_s(min_phys,max_phys,bcdef, &
     &           min_sub,max_sub,nsubs,subs_bc,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE define_subsbc_d(min_phys,max_phys,bcdef, &
     &           min_sub,max_sub,nsubs,subs_bc,info)
#endif
      !!! This routine defines the boundary conditions of the
      !!! subs on all processors. Thus subs that have
      !!! faces at the physical boundary are marked with +1 and
      !!! internal faces with the value 0. This information is
      !!! obtained by comparing the coordinates of the subvs with
      !!! the values of `min_phys` and `max_phys`.
      !!!
      !!! [WARNING]
      !!! If the user created the subs the comparison of floats in the present
      !!! routine might fail (round off errors) and the comparison should be
      !!! replaced by a `ABS(value-target) < epsilon` comparison.
      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_error
      USE ppm_module_alloc
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
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: min_phys
      !!! Min. extent of the physical domain
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: max_phys
      !!! Max. extent of the physical domain
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: bcdef
      !!! Boundary condition definition
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: min_sub
      !!! Min. extent of the sub domain
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: max_sub
      !!! Max. extent of the sub domain
      INTEGER                 , INTENT(IN   ) :: nsubs
      !!! The total number of (real) sub domains
      INTEGER , DIMENSION(:,:), POINTER       :: subs_bc
      !!! Boundary defintion of each sub
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER , DIMENSION(ppm_dim) :: ldc
      INTEGER                      :: i,iopt
      REAL(MK)                     :: t0
      REAL(MK)                     :: lmyeps
      !-------------------------------------------------------------------------
      !  Externals 
      !-------------------------------------------------------------------------
      
      !-------------------------------------------------------------------------
      !  Initialise 
      !-------------------------------------------------------------------------
      CALL substart('ppm_define_subs_bc',t0,info)
#if    __KIND == __SINGLE_PRECISION
      lmyeps = ppm_myepss
#elif  __KIND == __DOUBLE_PRECISION
      lmyeps = ppm_myepsd
#endif

      !-------------------------------------------------------------------------
      !  Allocate memory for subdomain bc flags
      !-------------------------------------------------------------------------
      iopt   = ppm_param_alloc_fit
      ldc(1) = 2*ppm_dim
      ldc(2) = MAX(nsubs,1)
      CALL ppm_alloc(subs_bc,ldc,iopt,info)
      IF (info.NE.0) THEN
          info = ppm_error_fatal
          CALL ppm_error(ppm_err_alloc,'ppm_define_subs_bc',     &
     &        'allocation of subs_bc failed',__LINE__,info)
          GOTO 9999
      ENDIF

      !-------------------------------------------------------------------------
      !  Loop over the global subs and compare their 
      !  coordinates with the physical boundary
      !-------------------------------------------------------------------------
      DO i=1,nsubs
         !----------------------------------------------------------------------
         !  compare the west boundary
         !----------------------------------------------------------------------
         IF (ABS(min_sub(1,i)-min_phys(1)) .LT. lmyeps*(max_sub(1,i)-min_sub(1,i))) THEN 
            subs_bc(1,i) = 1
         ELSE
            subs_bc(1,i) = 0
         ENDIF 

         !----------------------------------------------------------------------
         !  compare the east boundary
         !----------------------------------------------------------------------
         IF (ABS(max_sub(1,i)-max_phys(1)) .LT. lmyeps*(max_sub(1,i)-min_sub(1,i))) THEN
            subs_bc(2,i) = 1
         ELSE
            subs_bc(2,i) = 0
         ENDIF 

         !----------------------------------------------------------------------
         !  compare the south boundary
         !----------------------------------------------------------------------
         IF (ABS(min_sub(2,i)-min_phys(2)) .LT. lmyeps*(max_sub(2,i)-min_sub(2,i))) THEN
            subs_bc(3,i) = 1
         ELSE
            subs_bc(3,i) = 0
         ENDIF 

         !----------------------------------------------------------------------
         !  compare the north boundary
         !----------------------------------------------------------------------
         IF (ABS(max_sub(2,i)-max_phys(2)) .LT. lmyeps*(max_sub(2,i)-min_sub(2,i))) THEN
            subs_bc(4,i) = 1
         ELSE
            subs_bc(4,i) = 0
         ENDIF 

         !----------------------------------------------------------------------
         !  in three dimensions
         !----------------------------------------------------------------------
         IF (ppm_dim.EQ.3) THEN
            !-------------------------------------------------------------------
            !  compare the bottom boundary
            !-------------------------------------------------------------------
            IF (ABS(min_sub(3,i)-min_phys(3)) .LT. lmyeps*(max_sub(3,i)-min_sub(3,i))) THEN
               subs_bc(5,i) = 1
            ELSE
               subs_bc(5,i) = 0
            ENDIF 

            !-------------------------------------------------------------------
            !  compare the top boundary
            !-------------------------------------------------------------------
            IF (ABS(max_sub(3,i)-max_phys(3)) .LT. lmyeps*(max_sub(3,i)-min_sub(3,i))) THEN
               subs_bc(6,i) = 1
            ELSE
               subs_bc(6,i) = 0
            ENDIF 
         ENDIF 
      ENDDO

      !-------------------------------------------------------------------------
      !  Return 
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_define_subs_bc',t0,info)
      RETURN
#if   __KIND == __SINGLE_PRECISION
      END SUBROUTINE define_subsbc_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE define_subsbc_d
#endif
