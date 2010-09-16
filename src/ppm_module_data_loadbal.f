      !--*- f90 -*--------------------------------------------------------------
      !  Module       :              ppm_module_data_loadbal
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      MODULE ppm_module_data_loadbal
      !!! This module holds data used by the load balancing routines.
      !!!
      !!! [NOTE]
      !!! The variables declared in this modules should not be directly accessed
      !!! by the user. They are managed interally by the library.
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
