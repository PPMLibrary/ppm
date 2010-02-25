      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_define_subs_bc
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This routine defines the boundary conditions of the 
      !                 subs on all processors. Thus subs that have 
      !                 faces at the physical boundary are marked with +1 and
      !                 internal faces with the value 0. This information is 
      !                 obtained by comparing the coordinates of the subvs with
      !                 the values of min_phys and max_phys.
      !
      !  Input        : min_phys(:)  (F) : the min. extent of the physical domain
      !                 max_phys(:)  (F) : the max. extent of the physical domain
      !                 nsubs        (I) : the total number of (real) sub domains
      !                 bcdef(:)     (I) : boundary condition definition
      !               : min_sub(:,:) (F) : the min. extent of the sub domain
      !                 max_sub(:,:) (F) : the max. extent of the sub domain
      !
      !  Output       : subs_bc(:,:) (F) : boundary defintion of each sub
      !                 info         (I) : return status
      !
      !  Remarks      : WARNING: If the user created the subs the comparison
      !                 of floats in the present routine might fail (round off
      !                 errors) and the comparison should be replaced by a 
      !                 ABS(value-target).LT.epsilon comparison.
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_define_subs_bc.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.4  2006/09/05 08:01:26  pchatela
      !  Proper scaling for REAL comparisons
      !  Added module_alloc to ppm_decomp_boxsplit
      !
      !  Revision 1.3  2006/02/03 09:34:00  ivos
      !  Fixed bug 00015: ppm_subs_bc was only allocated and stored for the
      !  local subs in topo_store. Several mapping routines however need the
      !  info about all (global) subs.
      !  Changed subs_bc to hold now the GLOBAL subid and adjusted all
      !  occurrences.
      !
      !  Revision 1.2  2004/10/01 16:08:57  ivos
      !  Replaced REAL(ppm_kind_double) :: t0 with REAL(MK) t0.
      !
      !  Revision 1.1  2004/07/26 08:55:29  ivos
      !  Renamed.
      !
      !  Revision 1.3  2004/07/26 07:42:37  ivos
      !  Changed to single-interface modules and adapted all USE statements.
      !
      !  Revision 1.2  2004/04/07 11:20:16  ivos
      !  Corrected some copy-paste errors in comments and added lmyeps for
      !  REAL comparisons.
      !
      !  Revision 1.1  2004/04/05 12:23:04  walther
      !  First version.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KIND == __SINGLE_PRECISION
      SUBROUTINE ppm_define_subs_bc_s(min_phys,max_phys,bcdef, &
     &           min_sub,max_sub,nsubs,subs_bc,info)
#elif __KIND == __DOUBLE_PRECISION
      SUBROUTINE ppm_define_subs_bc_d(min_phys,max_phys,bcdef, &
     &           min_sub,max_sub,nsubs,subs_bc,info)
#endif
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
#include "ppm_define.h"
#ifdef __MPI
      INCLUDE 'mpif.h'
#endif
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK), DIMENSION(:)  , INTENT(IN   ) :: min_phys,max_phys
      INTEGER , DIMENSION(:  ), INTENT(IN   ) :: bcdef
      REAL(MK), DIMENSION(:,:), INTENT(IN   ) :: min_sub,max_sub
      INTEGER                 , INTENT(IN   ) :: nsubs
      INTEGER , DIMENSION(:,:), POINTER       :: subs_bc
      INTEGER                 , INTENT(  OUT) :: info
      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      INTEGER , DIMENSION(ppm_dim) :: ldc
      INTEGER               :: i,j,k,iopt
      REAL(MK)              :: t0
      REAL(MK)              :: lmyeps
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
      iopt   = ppm_param_alloc_grow_preserve
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
      END SUBROUTINE ppm_define_subs_bc_s
#elif __KIND == __DOUBLE_PRECISION
      END SUBROUTINE ppm_define_subs_bc_d
#endif
