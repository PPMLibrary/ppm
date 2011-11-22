  !----------------------------------------------------------
  ! symmetry ghosts in the x-direction.
  !----------------------------------------------------------
  
  !----------------------------------------------------------
  ! those at the west boundary.
  ! No '=' in the inequality checking.
  !----------------------------------------------------------
  
  xminf = min_phys(1) - ghostsize
  xmini = min_phys(1) 
  
  DO i = 1, Nall
     
     IF ( xp(1,i) > xminf .AND.xp(1,i) < xmini ) THEN
        xp(1,i)   = 2.0_MK*xminf - xp(1,i)
     END IF
  END DO
