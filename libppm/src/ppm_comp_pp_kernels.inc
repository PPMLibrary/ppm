      !-------------------------------------------------------------------------
      !  Include file for different PP kernels. This is included in all
      !  the pse_comp_pp* routines wherever kernel evaluations are needed.
      !  Only change this file to add/edit kernel functions.
      !
      !  INPUT:   dij is the squared distance between particles i and j
      !           All kernel parameters need to be passed in kpar(:).
      !           dx, dy and dz are the components of the inter particle vector. 
      !  OUTPUT:  eta is the kernel value for the given distance
      !
      !  Remark: Should SXF90 stop to unfold the loops (swap IF and DO)
      !  because the number of kernels grows too big, we will simply make
      !  several such .inc files and make several versions of the comp_pp
      !  routines with cpp if to choose the proper include file.
      !-------------------------------------------------------------------------
      !  $Log: ppm_comp_pp_kernels.inc,v $
      !  Revision 1.1.1.1  2007/07/13 10:18:54  ivos
      !  CBL version of the PPM library
      !
      !  Revision 1.4  2004/11/11 15:27:09  ivos
      !  Optimized PSE kernels for speed.
      !
      !  Revision 1.3  2004/10/13 16:42:11  davidch
      !  added support for 3d sph kernels
      !
      !  Revision 1.2  2004/07/29 15:56:48  hiebers
      !  added kernel_sph2d_p2 and kernel_dx_sph2d_p2
      !
      !  Revision 1.1  2004/07/23 12:57:31  ivos
      !  Initial implementation. Not tested.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   __KERNEL == __INTERNAL
      IF (kernel .EQ. ppm_param_kernel_fast3d) THEN
          !---------------------------------------------------------------------
          !  2nd order sph kernel 3d
          !     kpar(1)     cutoff^2
          !     kpar(2)     1/cutoff^9
          !     eta(delta^2) = 315/64/Pi/cutoff^9*(cutoff^2-r^2)^3
          !---------------------------------------------------------------------
          factor = kpar(1) - dij
          eta    = 1.5666814710608_MK * kpar(2)
          eta    = eta * factor * factor * factor
      ELSEIF (kernel .EQ. ppm_param_kernel_fast3d_dx) THEN
          !---------------------------------------------------------------------
          !  2nd order sph kernel 3d, dx
          !     kpar(1)     cutoff^2
          !     kpar(2)     1/cutoff^9
          !     eta(delta^2) = -945/32/Pi/cutoff^9*(cutoff^2-r^2)^2*x
          !---------------------------------------------------------------------
          factor = kpar(1) - dij
          eta    = -9.400088826365_MK * kpar(2) * dx
          eta    = eta * factor * factor
      ELSEIF (kernel .EQ. ppm_param_kernel_fast3d_dy) THEN
          !---------------------------------------------------------------------
          !  2nd order sph kernel 3d, dy
          !     kpar(1)     cutoff^2
          !     kpar(2)     1/cutoff^9
          !     eta(delta^2) = -945/32/Pi/cutoff^9*(cutoff^2-r^2)^2*y
          !---------------------------------------------------------------------
          factor = kpar(1) - dij
          eta    = -9.400088826365_MK * kpar(2) * dy
          eta    = eta * factor * factor
      ELSEIF (kernel .EQ. ppm_param_kernel_fast3d_dz) THEN
          !---------------------------------------------------------------------
          !  2nd order sph kernel 3d, dz
          !     kpar(1)     cutoff^2
          !     kpar(2)     1/cutoff^9
          !     eta(delta^2) = -945/32/Pi/cutoff^9*(cutoff^2-r^2)^2*z
          !---------------------------------------------------------------------
          factor = kpar(1) - dij
          eta    = -9.400088826365_MK * kpar(2) * dz
          eta    = eta * factor * factor
      ELSEIF (kernel .EQ. ppm_param_kernel_fast3d_lap) THEN
          !---------------------------------------------------------------------
          !  2nd order sph kernel 3d, laplacian
          !     kpar(1)     cutoff^2
          !     kpar(2)     1/cutoff^9
          ! eta(delta^2) = -945/32/Pi/cutoff^9*(cutoff^2-r^2)(3*cutoff^2-7*r^2)
          !---------------------------------------------------------------------
          factor  = kpar(1) - dij
          factor2 = 3.0_MK * kpar(1) - 7.0_MK * dij
          eta     = -9.400088826365_MK * kpar(2)
          eta     = eta * factor * factor2
      ELSEIF (kernel .EQ. ppm_param_kernel_laplace2d_p2) THEN
          !---------------------------------------------------------------------
          !  2nd order polynomial kernel 2D.
          !     kpar(1)     eps2inv  (1/eps**2)
          !     kpar(2)     correction*15.0*dh^2*eps^(-4)/pi^2
          !---------------------------------------------------------------------
          dij  = dij*kpar(1)
          dij2 = dij*dij
          dij4 = dij2*dij2
          dij5 = dij4*dij
          eta  = dij5 + 1.0_MK
          eta  = kpar(2)/eta
      ELSEIF (kernel .EQ. ppm_param_kernel_laplace3d_p2) THEN
          !---------------------------------------------------------------------
          !  2nd order polynomial kernel 3D.
          !     kpar(1)     eps2inv (1/eps**2)
          !     kpar(2)     correction*15.0*dh^3*eps^(-5)/pi^2
          !---------------------------------------------------------------------
          dij  = dij*kpar(1)
          dij2 = dij*dij
          dij4 = dij2*dij2
          dij5 = dij4*dij
          eta  = dij5 + 1.0_MK
          eta  = kpar(2)/eta
      ELSEIF (kernel .EQ. ppm_param_kernel_sph2d_p2 ) THEN
          !---------------------------------------------------------------------
          !  2nd order quartic spline kernel M5 2D.
          !     kpar(1)     eps (kernel width)
          !     kpar(2)     eps2 (eps**2)
          !     kpar(3)     eps2inv  (1/eps**2)
          !     kpar(4)     eps5inv  (1/eps**4)
          !     kpar(5)     dh3  (particle volume)
          !     kpar(6)     piinv  (1/PI, the math constant)
          !     kpar(7)     kappa (discretisation correction factor)
          !--------------------------------------------------------------------- 
          factor = kpar(3)*kpar(6)
          dij2    = dij*kpar(3)
          dij     = sqrt( dij2 )
  
          IF (dij.LT.2.5_MK) THEN
             IF (dij.LT.1.5_MK) THEN
            IF (dij.LT.0.5_MK) THEN
                   eta =  0.3_MK*dij2*dij2 - 0.75_MK*dij2 + 0.71875_MK
                ELSE
               eta = -0.2_MK*dij2*dij2 + dij2*dij - 1.5_MK*dij2 &
        &                                             + 0.25_MK*dij+0.6875_MK 
            ENDIF
         ELSE
        dij4  = 2.5_MK - dij
            eta   = 0.05_MK*dij4*dij4*dij4*dij4
             ENDIF
          ELSE 
                eta   = 0.0_MK
          ENDIF
          eta  = eta*factor*kpar(7)
      ELSEIF (kernel .EQ. ppm_param_kernel_dx_sph2d_p2 ) THEN
          !---------------------------------------------------------------------
          !  2nd order quartic spline kernel M5 2D for first Derivative d/dx.
          !     kpar(1)     eps (kernel width)
          !     kpar(2)     epsinv (1/eps)
          !     kpar(3)     eps2inv  (1/eps**2)
          !     kpar(4)     eps5inv  (1/eps**4)
          !     kpar(5)     dh3  (particle volume)
          !     kpar(6)     piinv  (1/PI, the math constant)
          !     kpar(7)     kappa (discretisation correction factor)
          !--------------------------------------------------------------------- 
          factor = kpar(3)*kpar(6)
          dij2    = dij*kpar(3)
          dij     = sqrt( dij2 )
  
          IF ( dij.LT.2.5_MK .AND. dij>0.0_MK ) THEN
             IF (dij.LT.1.5_MK) THEN
                IF (dij.LT.0.5_MK) THEN
                   eta =  1.2_MK*dij2*dij - 1.5_MK*dij
                ELSE
                   eta = -0.8_MK*dij2*dij + 3_MK*dij2 - 3_MK*dij + 0.25_MK 
                ENDIF
             ELSE
                dij4  =  2.5_MK - dij
        eta   = -0.2_MK*dij4*dij4*dij4
         ENDIF
          ELSE 
        eta   = 0.0_MK
      ENDIF
         eta  = eta*factor*kpar(7)*kpar(2)/dij*dx
!          eta  = eta*factor*kpar(7)*kpar(2)/dij
       ELSE
          !---------------------------------------------------------------------
          !  This ELSE is needed to avoid funny action when an invalid
          !  kernel is specified (and to make SXF90 unfold the loops).
          !---------------------------------------------------------------------
          eta = 0.0_MK
      ENDIF
#elif __KERNEL == __USER_FUNCTION
      eta = kernel(dij,kpar)
#elif __KERNEL == __LOOKUP_TABLE
      idx = INT(dij*kpar)
      eta = kernel(idx)
#endif
