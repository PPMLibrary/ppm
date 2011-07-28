      !-------------------------------------------------------------------------
      !  Module       :              ppm_module_data_loadbal
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for the routines
      !                 concerned with load balancing. They are callable by
      !                 the user.
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_data_loadbal.f,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:57  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.2  2004/07/26 11:48:10  ivos
      !  Fixes to make it compile.
      !
      !  Revision 1.1  2004/07/26 07:28:17  ivos
      !  Initial implementation. Originated from splitting the old ppm
      !  modules.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      MODULE ppm_module_data_loadbal

         !----------------------------------------------------------------------
         !  Modules
         !----------------------------------------------------------------------
         USE ppm_module_data, ONLY: ppm_kind_double
         PRIVATE :: ppm_kind_double

         !----------------------------------------------------------------------
         !  Timing and load statistics
         !----------------------------------------------------------------------
         ! estimated cost of redefining the topology
         REAL(ppm_kind_double) :: ppm_loadbal_decomp_cost = 0.0_ppm_kind_double
         ! counter of how many times decomp cost was measured (N)
         INTEGER               :: ppm_loadbal_dcn = 0
         ! running sum of (maxtime - avgtime)
         REAL(ppm_kind_double) :: ppm_loadbal_runsum = 0.0_ppm_kind_double
         ! old value of the SAR function
         REAL(ppm_kind_double) :: ppm_loadbal_old_sar = -1.0_ppm_kind_double

      END MODULE ppm_module_data_loadbal
