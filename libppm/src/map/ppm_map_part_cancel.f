      !-------------------------------------------------------------------------
      !  Subroutine   :                 ppm_map_part_cancel
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  ETH Zurich
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      SUBROUTINE ppm_map_part_cancel(info)
      !!! This routine cancels an ongoing mapping calling sequence 
      !!! and resets the internal mapping state variables.
      !!!
      !!! [NOTE]
      !!! one needs to reset the arrays since thay are alloc_fit-ed the next time a
      !!! mapping is called. Their deallocation only happens in
      !!! ppm_finalize anyway.
      ! TODO: Fix this awful docstring!

      !-------------------------------------------------------------------------
      !  Modules
      !-------------------------------------------------------------------------
      USE ppm_module_substart
      USE ppm_module_substop
      USE ppm_module_data
      IMPLICIT NONE


      !-------------------------------------------------------------------------
      !  Arguments
      !-------------------------------------------------------------------------
      INTEGER                 , INTENT(  OUT) :: info
      !!! Returns status, 0 upon success
      !-------------------------------------------------------------------------
      !  Local Variables
      !-------------------------------------------------------------------------
      REAL(ppm_kind_double)                :: t0
      !-------------------------------------------------------------------------
      !  Initialise
      !-------------------------------------------------------------------------
      CALL substart('ppm_map_part_cancel',t0,info)



      ppm_buffer_set    = 0
      ppm_nsendbuffer   = 0
      ppm_nrecvbuffer   = 0

      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
 9999 CONTINUE
      CALL substop('ppm_map_part_cancel',t0,info)
      RETURN
      END SUBROUTINE ppm_map_part_cancel
