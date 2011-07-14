  !----------------------------------------------------------
  ! (Re)allocate memory for the symmetry ghosts.
  !----------------------------------------------------------
  
  iopt   = ppm_param_alloc_grow_preserve
  ldu(1) = ppm_dim
  ldu(2) = 2*nghostplus
  CALL ppm_alloc(xt,ldu,iopt,info) 
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
  
  k  = nghostplus

  !----------------------------------------------------------
  ! those at the west boundary.
  ! No '=' in the inequality checking.
  !----------------------------------------------------------
  
  xminf = min_phys(1) 
  xmini = min_phys(1) + ghostsize
  
  DO i = 1, nghostplus
     
     IF ( xt(1,i) > xminf .AND.xt(1,i) < xmini ) THEN
        k         = k + 1
        ighost(k) = ighost(i)
        xt(1,k)   = 2.0_MK*xminf - xt(1,i)
        xt(2,k)   = xt(2,i)
        IF ( ppm_dim == 3 ) THEN
           xt(3,k)   = xt(3,i)
        END IF
     END IF
  END DO
  
  !----------------------------------------------------------
  ! update the ghost counter.
  !-----------------------------------------------------------
  
  nghostplus = k
  
