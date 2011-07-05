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
  ! Symmetry ghosts in the z-direction.
  !----------------------------------------------------------

  k = nghostplus

  !----------------------------------------------------------
  ! those at upper boundary, 
  ! even using symmetry inter-communication.
  ! No '=' in the inequality checking.
  !----------------------------------------------------------
  
  zmaxf = max_phys(3) 
  zmaxi = max_phys(3) - ghostsize
  
  DO i = 1, nghostplus
     
     IF ( xt(3,i) > zmaxi .AND. &
          xt(3,i) < zmaxf ) THEN
        
        k         = k + 1
        ighost(k) = ighost(i)
        xt(1,k)   = xt(1,i)
        xt(2,k)   = xt(2,i) 
        xt(3,k)   = 2.0_MK*zmaxf - xt(3,i)
        
     END IF
     
  END DO
  
  !----------------------------------------------------------
  ! update the ghost counter.
  !----------------------------------------------------------
  
  nghostplus = k
  
  
