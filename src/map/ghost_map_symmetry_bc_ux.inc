  !----------------------------------------------------------
  ! (Re)allocate memory for the symmetry ghosts.
  !----------------------------------------------------------
  
  iopt   = ppm_param_alloc_grow_preserve
  ldu(1) = ppm_dim
  ldu(2) = 2*nghostplus
  CALL ppm_alloc(xt,ldu,iopt,info) 
  CALL ppm_alloc(xt_offset,ldu,iopt,info) 
  IF (info.NE.0) THEN
     info = ppm_error_fatal
     CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get', &
          'xt',__LINE__,info)
     GOTO 9999
  ENDIF
  
  ldu(1) = ldu(2)
  CALL ppm_alloc(ighost,ldu,iopt,info) 
  IF (info.NE.0) THEN
     info = ppm_error_fatal
     CALL ppm_error(ppm_err_alloc,'ppm_map_part_ghost_get', &
          'ighost',__LINE__,info)
     GOTO 9999
  ENDIF
  
  !----------------------------------------------------------
  ! symmetry ghosts in the x-direction.
  !----------------------------------------------------------

  k = nghostplus
  
  !----------------------------------------------------------
  ! the east bc, even we are using symmetry.
  ! No '=' in the inequality checking.
  !----------------------------------------------------------
  
  xmaxf = max_phys(1)
#if    __MODE == __HOM 
  xmaxi = max_phys(1) - ghostsize
#elif  __MODE == __INHOM
   !----------------------------------------------------------------
   !  For inhomogenous ghostlayers, we use the entire max
   !----------------------------------------------------------------
#if    __KIND == __SINGLE_PRECISION  
  xmaxi = max_phys(1) - topo%ghostsizes(1)
#else
  xmaxi = max_phys(1) - topo%ghostsized(1)
#endif
#endif
  
  
  DO i = 1, nghostplus
     
     IF ( xt(1,i) > xmaxi .AND. xt(1,i) < xmaxf) THEN
        k         = k + 1
        ighost(k) = ighost(i)
        xt(1,k)   = 2.0_MK*xmaxf - xt(1,i)
        xt(2,k)   = xt(2,i)
        xt_offset(1,k) = 2.0_MK*(xmaxf-xt(1,i))
        xt_offset(2,k) = xt_offset(2,i)
        IF ( ppm_dim == 3) THEN
           xt(3,k)   = xt(3,i)
        END IF
     END IF
  END DO
  
  !----------------------------------------------------------
  ! update the ghost counter
  !----------------------------------------------------------
  
  nghostplus = k
  
  
